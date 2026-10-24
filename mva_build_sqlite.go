// Parallel BAM→SQLite builder: reads gene-tagged BAM/SAM, extracts (cb, umi, gene),
// and writes to SQLite with the same schema as build_sqlite4.py.
//
// Design:
//   - BAM input: pipes through `samtools view` (C-level BAM decoding)
//   - SAM/SAM.gz input: direct file reading
//   - Pipeline: reader goroutine → parser pool (N workers) → SQLite writer (1 goroutine)
//   - Single transaction for all inserts (10-100x faster than per-batch commit)
//   - Stats via SQL COUNT queries at end
//
// Speedup vs Python:
//   - Parsing: Go strings.Split + parallel goroutines, no GIL
//   - SQLite: single transaction (vs per-batch commit in Python = multiple fsync)
//   - I/O: pipelined reader/parser/writer
//   - Expected: 3-8x faster overall
//
// Build:
//   go build -o build_sqlite build_sqlite4.go
//
// Usage:
//   build_sqlite -i gene.sorted.bam -o output/gene.molecules.sqlite [options]
//   build_sqlite -i gene.sam --qname-simple -o output.sqlite
//   build_sqlite -i gene.sorted.bam --use-cb-tag -o output.sqlite -t 8
//   build_sqlite --merge db1.sqlite,db2.sqlite -o merged.sqlite
//
// Requires: samtools in PATH (for BAM input). No CGO/gcc needed (pure Go SQLite).

package main

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"time"

	_ "modernc.org/sqlite"
)

// ─── Types ───────────────────────────────────────────────────────────────────

// Molecule is a single (cb, umi, gene) record.
type Molecule struct {
	CB, UMI, Gene string
}

// WorkerStats tracks per-worker parsing statistics.
type WorkerStats struct {
	Total, BadSAM, Unmapped, Secondary, LowMapq, NoGene, NoCbUmi, Counted int64
}

func (s *WorkerStats) Add(o WorkerStats) {
	s.Total += o.Total
	s.BadSAM += o.BadSAM
	s.Unmapped += o.Unmapped
	s.Secondary += o.Secondary
	s.LowMapq += o.LowMapq
	s.NoGene += o.NoGene
	s.NoCbUmi += o.NoCbUmi
	s.Counted += o.Counted
}

// ParsedChunk is a batch of parsed molecules + stats.
type ParsedChunk struct {
	Molecules []Molecule
	Stats     WorkerStats
}

// ─── Input opening ───────────────────────────────────────────────────────────

type readCloser struct {
	io.Reader
	closeFn func() error
}

func (rc *readCloser) Close() error { return rc.closeFn() }

type commandReader struct {
	io.Reader
	cmd *exec.Cmd
}

func (c *commandReader) Close() error {
	if closer, ok := c.Reader.(io.Closer); ok {
		closer.Close()
	}
	return c.cmd.Wait()
}

func openInput(path string) (io.ReadCloser, error) {
	if path == "-" {
		return io.NopCloser(os.Stdin), nil
	}

	lower := strings.ToLower(path)
	if strings.HasSuffix(lower, ".bam") || strings.HasSuffix(lower, ".cram") {
		// BAM/CRAM: pipe through samtools for C-level decoding
		cmd := exec.Command("samtools", "view", path)
		stdout, err := cmd.StdoutPipe()
		if err != nil {
			return nil, err
		}
		if err := cmd.Start(); err != nil {
			return nil, fmt.Errorf("samtools not found or failed: %w", err)
		}
		return &commandReader{stdout, cmd}, nil
	}

	// SAM (possibly gzip-compressed)
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	buf := make([]byte, 2)
	n, _ := f.Read(buf)
	f.Seek(0, 0)
	if n == 2 && buf[0] == 0x1f && buf[1] == 0x8b {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		return &readCloser{gz, func() error { gz.Close(); return f.Close() }}, nil
	}
	return f, nil
}

// ─── CB/UMI extraction ───────────────────────────────────────────────────────

type cbUmiExtractor struct {
	re     *regexp.Regexp
	cbIdx  int
	umiIdx int
	useTag bool
}

func newExtractor(qnameRegex string, useCbTag, qnameSimple bool) *cbUmiExtractor {
	if useCbTag {
		return &cbUmiExtractor{useTag: true}
	}

	var pattern string
	if qnameSimple {
		pattern = `_([^_]+)_(.+)$`
	} else if qnameRegex != "" {
		pattern = qnameRegex
	} else {
		pattern = `_(?P<cb>[ACGTN]+)_(?P<umi>[ACGTN]+)$`
	}

	re := regexp.MustCompile(pattern)
	names := re.SubexpNames()
	cbIdx, umiIdx := -1, -1
	for i, name := range names {
		switch name {
		case "cb":
			cbIdx = i
		case "umi":
			umiIdx = i
		}
	}
	// Fallback: if no named groups, assume group 1=cb, group 2=umi
	if cbIdx < 0 {
		cbIdx = 1
	}
	if umiIdx < 0 {
		umiIdx = 2
	}
	return &cbUmiExtractor{re: re, cbIdx: cbIdx, umiIdx: umiIdx}
}

func (e *cbUmiExtractor) extract(qname string, fields []string) (cb, umi string, ok bool) {
	if e.useTag {
		for i := 11; i < len(fields); i++ {
			f := fields[i]
			if strings.HasPrefix(f, "CB:Z:") {
				cb = f[5:]
			} else if strings.HasPrefix(f, "UB:Z:") || strings.HasPrefix(f, "UR:Z:") {
				umi = f[5:]
			}
		}
		return cb, umi, cb != "" && umi != ""
	}

	m := e.re.FindStringSubmatch(qname)
	if m == nil {
		return "", "", false
	}
	return m[e.cbIdx], m[e.umiIdx], true
}

// ─── Gene extraction ─────────────────────────────────────────────────────────

func geneFromRname(rname, delimiter string) string {
	if rname == "*" || !strings.Contains(rname, delimiter) {
		return ""
	}
	idx := strings.LastIndex(rname, delimiter)
	return rname[idx+1:]
}

// ─── SAM line parsing ────────────────────────────────────────────────────────

func parseLines(lines []string, ext *cbUmiExtractor, geneDelimiter string, minMapq int, skipSecondary bool) ([]Molecule, WorkerStats) {
	var mols []Molecule
	var stats WorkerStats

	for _, line := range lines {
		stats.Total++

		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			stats.BadSAM++
			continue
		}

		flag, err1 := strconv.Atoi(fields[1])
		mapq, err2 := strconv.Atoi(fields[4])
		if err1 != nil || err2 != nil {
			stats.BadSAM++
			continue
		}

		if flag&0x4 != 0 {
			stats.Unmapped++
			continue
		}
		if skipSecondary && (flag&0x100 != 0 || flag&0x800 != 0) {
			stats.Secondary++
			continue
		}
		if mapq < minMapq {
			stats.LowMapq++
			continue
		}

		gene := geneFromRname(fields[2], geneDelimiter)
		if gene == "" {
			stats.NoGene++
			continue
		}

		cb, umi, ok := ext.extract(fields[0], fields)
		if !ok {
			stats.NoCbUmi++
			continue
		}

		mols = append(mols, Molecule{cb, umi, gene})
		stats.Counted++
	}

	return mols, stats
}

// ─── SQLite building ─────────────────────────────────────────────────────────

func buildSQLite(
	inputPath string,
	ext *cbUmiExtractor,
	geneDelimiter string,
	minMapq int,
	skipSecondary bool,
	sqlitePath string,
	batchSize int,
	progressEvery int,
	threads int,
) error {
	// Ensure output directory exists
	os.MkdirAll(filepath.Dir(sqlitePath), 0755)

	// Open SQLite with single connection (WAL mode)
	db, err := sql.Open("sqlite", sqlitePath)
	if err != nil {
		return err
	}
	defer db.Close()
	db.SetMaxOpenConns(1)

	db.Exec("PRAGMA journal_mode=WAL")
	db.Exec("PRAGMA synchronous=OFF")
	db.Exec("PRAGMA temp_store=MEMORY")
	db.Exec("PRAGMA cache_size=-200000")
	db.Exec("DROP TABLE IF EXISTS molecules")
	db.Exec("CREATE TABLE molecules (cb TEXT NOT NULL, umi TEXT NOT NULL, gene TEXT NOT NULL)")

	// Open input
	rc, err := openInput(inputPath)
	if err != nil {
		return err
	}
	defer rc.Close()

	reader := bufio.NewReaderSize(rc, 4*1024*1024)

	// Pipeline: reader → chunkCh → parser pool → resultCh → writer
	chunkSize := 50000
	if threads <= 1 {
		chunkSize = 100000
	}

	chunkCh := make(chan []string, threads*2)
	resultCh := make(chan ParsedChunk, threads*2)

	// ─── Reader goroutine: reads lines, chunks them ───────────────────────
	go func() {
		var lines []string
		for {
			line, err := reader.ReadString('\n')
			if len(line) > 0 {
				line = strings.TrimRight(line, "\r\n")
				if line != "" && !strings.HasPrefix(line, "@") {
					lines = append(lines, line)
					if len(lines) >= chunkSize {
						chunkCh <- lines
						lines = nil
					}
				}
			}
			if err != nil {
				break
			}
		}
		if len(lines) > 0 {
			chunkCh <- lines
		}
		close(chunkCh)
	}()

	// ─── Parser workers: parse chunks in parallel ──────────────────────────
	var wg sync.WaitGroup
	for i := 0; i < threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for lines := range chunkCh {
				mols, stats := parseLines(lines, ext, geneDelimiter, minMapq, skipSecondary)
				resultCh <- ParsedChunk{mols, stats}
			}
		}()
	}

	// Close resultCh after all parsers finish
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// ─── Writer: single transaction, batch inserts ────────────────────────
	started := time.Now()
	var totalStats WorkerStats

	tx, err := db.Begin()
	if err != nil {
		return err
	}
	stmt, err := tx.Prepare("INSERT INTO molecules VALUES (?, ?, ?)")
	if err != nil {
		tx.Rollback()
		return err
	}

	batch := make([]Molecule, 0, batchSize)
	lastProgressBatch := int64(0)

	for pc := range resultCh {
		totalStats.Add(pc.Stats)
		batch = append(batch, pc.Molecules...)

		// Flush full batches
		for len(batch) >= batchSize {
			for _, m := range batch[:batchSize] {
				if _, err := stmt.Exec(m.CB, m.UMI, m.Gene); err != nil {
					stmt.Close()
					tx.Rollback()
					return fmt.Errorf("SQLite insert error: %w", err)
				}
			}
			batch = batch[batchSize:]
		}

		// Progress: print roughly every progressEvery reads
		if progressEvery > 0 {
			currentBatch := totalStats.Total / int64(progressEvery)
			if currentBatch > lastProgressBatch {
				lastProgressBatch = currentBatch
				elapsed := time.Since(started).Seconds()
				fmt.Fprintf(os.Stderr,
					"\r[sqlite] reads=%d counted=%d rate=%.0f/s",
					totalStats.Total, totalStats.Counted,
					float64(totalStats.Total)/elapsed)
			}
		}
	}

	// Insert remaining
	for _, m := range batch {
		stmt.Exec(m.CB, m.UMI, m.Gene)
	}
	stmt.Close()
	if err := tx.Commit(); err != nil {
		return err
	}

	// Create indexes after all inserts (faster than incremental)
	fmt.Fprintf(os.Stderr, "\n[sqlite] creating indexes...\n")
	db.Exec("CREATE INDEX idx_molecules_cb ON molecules(cb)")
	db.Exec("CREATE INDEX idx_molecules_cb_gene ON molecules(cb, gene)")

	// Query stats
	var nMolecules, nBarcodes, nGenes int
	db.QueryRow("SELECT COUNT(*) FROM molecules").Scan(&nMolecules)
	db.QueryRow("SELECT COUNT(DISTINCT cb) FROM molecules").Scan(&nBarcodes)
	db.QueryRow("SELECT COUNT(DISTINCT gene) FROM molecules").Scan(&nGenes)

	elapsed := time.Since(started).Seconds()
	fmt.Fprintf(os.Stderr,
		"\n[sqlite] done reads=%d counted=%d molecules=%d barcodes=%d genes=%d time=%.1fs\n",
		totalStats.Total, totalStats.Counted, nMolecules, nBarcodes, nGenes, elapsed)

	// Print summary to stdout (matches Python output format)
	fmt.Printf("output: %s\n", sqlitePath)
	fmt.Printf("barcodes: %d\n", nBarcodes)
	fmt.Printf("molecules: %d\n", nMolecules)
	fmt.Printf("genes: %d\n", nGenes)
	fmt.Println("--- read filtering breakdown ---")
	fmt.Printf("  reads_total: %d\n", totalStats.Total)
	fmt.Printf("  reads_bad_sam: %d\n", totalStats.BadSAM)
	fmt.Printf("  reads_unmapped: %d\n", totalStats.Unmapped)
	fmt.Printf("  reads_secondary_or_supplementary: %d\n", totalStats.Secondary)
	fmt.Printf("  reads_low_mapq: %d\n", totalStats.LowMapq)
	fmt.Printf("  reads_no_gene: %d\n", totalStats.NoGene)
	fmt.Printf("  reads_no_cb_umi: %d\n", totalStats.NoCbUmi)
	fmt.Printf("  reads_counted: %d\n", totalStats.Counted)

	return nil
}

// ─── Merge ───────────────────────────────────────────────────────────────────

func mergeSQLite(inputs []string, output string) error {
	if _, err := os.Stat(output); err == nil {
		os.Remove(output)
	}
	os.MkdirAll(filepath.Dir(output), 0755)

	db, err := sql.Open("sqlite", output)
	if err != nil {
		return err
	}
	defer db.Close()
	db.SetMaxOpenConns(1)

	db.Exec("PRAGMA journal_mode=WAL")
	db.Exec("PRAGMA synchronous=OFF")
	db.Exec("PRAGMA temp_store=MEMORY")
	db.Exec("PRAGMA cache_size=-200000")
	db.Exec("CREATE TABLE molecules (cb TEXT NOT NULL, umi TEXT NOT NULL, gene TEXT NOT NULL)")

	started := time.Now()
	total := 0

	for _, path := range inputs {
		absPath, _ := filepath.Abs(path)
		// Escape single quotes in path
		escapedPath := strings.ReplaceAll(absPath, "'", "''")
		_, err := db.Exec(fmt.Sprintf("ATTACH DATABASE 'file:%s?mode=ro' AS src_db", escapedPath))
		if err != nil {
			return fmt.Errorf("ATTACH failed for %s: %w", path, err)
		}
		_, err = db.Exec("INSERT INTO molecules SELECT cb, umi, gene FROM src_db.molecules")
		if err != nil {
			return fmt.Errorf("INSERT FROM %s: %w", path, err)
		}
		var count int
		db.QueryRow("SELECT COUNT(*) FROM src_db.molecules").Scan(&count)
		db.Exec("DETACH DATABASE src_db")
		total += count
		fmt.Fprintf(os.Stderr, "[merge] %s: %d molecules\n", path, count)
	}

	db.Exec("CREATE INDEX idx_molecules_cb ON molecules(cb)")
	db.Exec("CREATE INDEX idx_molecules_cb_gene ON molecules(cb, gene)")

	elapsed := time.Since(started).Seconds()
	fmt.Fprintf(os.Stderr, "[merge] done: %d inputs, %d total molecules, time=%.1fs\n",
		len(inputs), total, elapsed)
	fmt.Printf("output: %s\n", output)
	fmt.Printf("molecules: %d\n", total)

	return nil
}

// ─── Main ─────────────────────────────────────────────────────────────────────

func main() {
	var (
		inputPath     string
		outputPath    string
		qnameRegex    string
		useCbTag      bool
		qnameSimple   bool
		geneDelimiter string
		minMapq       int
		keepSecondary bool
		batchSize     int
		progressEvery int
		threads       int
		mergeStr      string
	)

	flag.StringVar(&inputPath, "i", "", "gene-tagged SAM/BAM (BAM requires samtools in PATH)")
	flag.StringVar(&outputPath, "o", "", "output SQLite path [required]")
	flag.StringVar(&qnameRegex, "qname-regex", `_(?P<cb>[ACGTN]+)_(?P<umi>[ACGTN]+)$`,
		"regex with named groups cb and umi")
	flag.BoolVar(&useCbTag, "use-cb-tag", false, "use CB:Z: and UB:Z:/UR:Z: SAM tags instead of qname regex")
	flag.BoolVar(&qnameSimple, "qname-simple", false, "simple regex _([^_]+)_(.+)$ (overrides --qname-regex)")
	flag.StringVar(&geneDelimiter, "gene-delimiter", "#", "delimiter between contig and gene in RNAME")
	flag.IntVar(&minMapq, "min-mapq", 0, "minimum MAPQ to count")
	flag.BoolVar(&keepSecondary, "keep-secondary", false, "also count secondary/supplementary alignments")
	flag.IntVar(&batchSize, "batch-size", 200000, "SQLite insert batch size")
	flag.IntVar(&progressEvery, "progress-every", 1000000, "print progress every N reads (0 disables)")
	flag.IntVar(&threads, "t", 4, "parser goroutines")
	flag.StringVar(&mergeStr, "merge", "", "comma-separated SQLite DBs to merge into -o output")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Parallel BAM→SQLite builder (Go version of build_sqlite4.py)\n\n")
		fmt.Fprintf(os.Stderr, "Usage:\n")
		fmt.Fprintf(os.Stderr, "  build_sqlite -i gene.sorted.bam -o output.sqlite [options]\n")
		fmt.Fprintf(os.Stderr, "  build_sqlite -i gene.sam --qname-simple -o output.sqlite\n")
		fmt.Fprintf(os.Stderr, "  build_sqlite -i gene.sorted.bam --use-cb-tag -o output.sqlite -t 8\n")
		fmt.Fprintf(os.Stderr, "  build_sqlite --merge db1.sqlite,db2.sqlite -o merged.sqlite\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	// Merge mode
	if mergeStr != "" {
		if outputPath == "" {
			fmt.Fprintln(os.Stderr, "[ERROR] -o is required for --merge")
			os.Exit(1)
		}
		inputs := strings.Split(mergeStr, ",")
		for i := range inputs {
			inputs[i] = strings.TrimSpace(inputs[i])
		}
		if err := mergeSQLite(inputs, outputPath); err != nil {
			fmt.Fprintf(os.Stderr, "[ERROR] %v\n", err)
			os.Exit(1)
		}
		return
	}

	// Build mode
	if inputPath == "" || outputPath == "" {
		fmt.Fprintln(os.Stderr, "[ERROR] -i and -o are required")
		flag.Usage()
		os.Exit(1)
	}
	if batchSize < 1000 {
		batchSize = 200000
	}
	if threads < 1 {
		threads = 1
	}

	// Build extractor
	ext := newExtractor(qnameRegex, useCbTag, qnameSimple)

	fmt.Fprintf(os.Stderr, "[INFO] input=%s output=%s threads=%d cb_tag=%v qname_simple=%v\n",
		inputPath, outputPath, threads, useCbTag, qnameSimple)

	if err := buildSQLite(inputPath, ext, geneDelimiter, minMapq, !keepSecondary,
		outputPath, batchSize, progressEvery, threads); err != nil {
		fmt.Fprintf(os.Stderr, "\n[ERROR] %v\n", err)
		os.Exit(1)
	}
}
