// Parallel UMI+Barcode extraction: extract BC/UMI from R1 and append to R2 header.
//
// Design:
//   - Reader goroutine reads R1+R2 in large chunks, sends to worker pool
//   - Worker goroutines do BC/UMI extraction and format new headers
//   - Results are reassembled in order and written to output
//   - Supports both plain FASTQ and gzip (.fq.gz)
//   - Verifies read ID (/1 vs /2) match for each pair; aborts on mismatch
//
// Output header format (umi_tools extract style):
//   @READ_ID_BC_UMI
//   e.g. @LH29866:259:22L7IQLT4:8:1108:7974:29714_TCATCGTATCGGTCCCTAGT_ATAGCGAA
//
// Usage:
//   umi_extract -1 r1.fq.gz -2 r2.fq.gz -o r2_extracted.fq.gz \
//       --bc-start 1 --bc-len 20 --umi-start 21 --umi-len 8 -t 12
//
// Default: bc=1-20bp, umi=21-28bp, 12 threads, chunk=200000 records

package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"
	"sync"
	"time"
)

// ─── Data types ──────────────────────────────────────────────────────────────

// Record holds a paired FASTQ record (minimal fields needed).
type Record struct {
	H1, S1 string // R1 header (with @) and sequence
	H2, S2 string // R2 header (with @) and sequence
	P2, Q2 string // R2 "+" line and quality
}

// OutRecord is a processed record ready for output.
type OutRecord struct {
	Header, Seq, Plus, Qual string
}

// Chunk is a batch of records with a sequential index for ordered output.
type Chunk struct {
	Index   int
	Records []Record
}

// Result is a processed chunk.
type Result struct {
	Index   int
	Records []OutRecord
	Err     error
}

// ─── File I/O helpers ────────────────────────────────────────────────────────

// readCloser wraps a gzip.Reader + os.File so Close() closes both.
type readCloser struct {
	io.Reader
	closeFn func() error
}

func (rc *readCloser) Close() error { return rc.closeFn() }

func openReader(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		return &readCloser{gz, func() error { gz.Close(); return f.Close() }}, nil
	}
	return f, nil
}

// openWriter returns a buffered writer and a cleanup function.
// If path ends with .gz, output is gzip-compressed (level 4, same as Python).
func openWriter(path string) (*bufio.Writer, func(), error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, nil, err
	}
	var w io.Writer = f
	var gz *gzip.Writer
	if strings.HasSuffix(path, ".gz") {
		gz, err = gzip.NewWriterLevel(f, 4)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		w = gz
	}
	bw := bufio.NewWriterSize(w, 4*1024*1024)
	cleanup := func() {
		bw.Flush()
		if gz != nil {
			gz.Close()
		}
		f.Close()
	}
	return bw, cleanup, nil
}

// ─── FASTQ reading ───────────────────────────────────────────────────────────

// readLine reads one line, trimming trailing \r\n. Returns (line, cleanEOF).
func readLine(r *bufio.Reader) (string, bool, error) {
	line, err := r.ReadString('\n')
	if err == io.EOF {
		trimmed := strings.TrimRight(line, "\r\n")
		if len(trimmed) == 0 {
			return "", true, nil // clean EOF
		}
		return trimmed, false, nil // last line without newline
	}
	if err != nil {
		return "", false, err
	}
	return strings.TrimRight(line, "\r\n"), false, nil
}

// readRequiredLine reads one line; returns error on unexpected EOF.
func readRequiredLine(r *bufio.Reader) (string, error) {
	line, eof, err := readLine(r)
	if err != nil {
		return "", err
	}
	if eof {
		return "", io.ErrUnexpectedEOF
	}
	return line, nil
}

// readChunk reads up to maxRecords paired FASTQ records from r1 and r2.
// Returns (records, cleanEOF, error).
func readChunk(r1, r2 *bufio.Reader, maxRecords int) ([]Record, bool, error) {
	records := make([]Record, 0, maxRecords)
	for i := 0; i < maxRecords; i++ {
		h1, eof, err := readLine(r1)
		if err != nil {
			return nil, false, err
		}
		if eof {
			return records, true, nil
		}

		s1, err := readRequiredLine(r1)
		if err != nil {
			return nil, false, fmt.Errorf("R1 unexpected EOF reading seq: %v", err)
		}
		readRequiredLine(r1) // +
		readRequiredLine(r1) // qual

		h2, eof2, err := readLine(r2)
		if err != nil {
			return nil, false, err
		}
		if eof2 {
			return nil, false, fmt.Errorf("R2 has fewer records than R1")
		}
		s2, err := readRequiredLine(r2)
		if err != nil {
			return nil, false, fmt.Errorf("R2 unexpected EOF reading seq: %v", err)
		}
		p2, _ := readRequiredLine(r2)
		q2, _ := readRequiredLine(r2)

		records = append(records, Record{h1, s1, h2, s2, p2, q2})
	}
	return records, false, nil
}

// ─── Chunk processing ────────────────────────────────────────────────────────

// processChunk extracts BC+UMI from R1 seq and builds new R2 headers.
func processChunk(chunk Chunk, bcStart, bcLen, umiStart, umiLen int) Result {
	out := make([]OutRecord, 0, len(chunk.Records))
	bcBegin := bcStart - 1   // 1-based → 0-based
	umiBegin := umiStart - 1 // 1-based → 0-based

	for _, rec := range chunk.Records {
		// Strip @ prefix, keep only the part before first space
		id1 := strings.TrimLeft(rec.H1, "@")
		if idx := strings.IndexByte(id1, ' '); idx > 0 {
			id1 = id1[:idx]
		}

		id2 := strings.TrimLeft(rec.H2, "@")
		if idx := strings.IndexByte(id2, ' '); idx > 0 {
			id2 = id2[:idx]
		}

		if id1 != id2 {
			return Result{
				Index: chunk.Index,
				Err:   fmt.Errorf("read ID mismatch: '%s' != '%s'", id1, id2),
			}
		}

		// Extract BC and UMI from R1 seq (1-based → 0-based slice)
		bc := safeSlice(rec.S1, bcBegin, bcLen)
		umi := safeSlice(rec.S1, umiBegin, umiLen)

		out = append(out, OutRecord{
			Header: "@" + id2 + "_" + bc + "_" + umi,
			Seq:    rec.S2,
			Plus:   rec.P2,
			Qual:   rec.Q2,
		})
	}

	return Result{Index: chunk.Index, Records: out}
}

// safeSlice returns s[start:start+len], clamped to string bounds.
func safeSlice(s string, start, length int) string {
	if start >= len(s) {
		return ""
	}
	end := start + length
	if end > len(s) {
		end = len(s)
	}
	return s[start:end]
}

// ─── Main ────────────────────────────────────────────────────────────────────

func defaultThreads() int {
	n := runtime.NumCPU()
	if n > 12 {
		n = 12
	}
	return n
}

func main() {
	r1Path := flag.String("1", "", "R1 FASTQ (.gz supported) [required]")
	r2Path := flag.String("2", "", "R2 FASTQ (.gz supported) [required]")
	outPath := flag.String("o", "", "Output R2 with UMI in header [required]")
	bcStart := flag.Int("bc-start", 1, "1-based BC start (default: 1)")
	bcLen := flag.Int("bc-len", 20, "BC length (default: 20)")
	umiStart := flag.Int("umi-start", 21, "1-based UMI start (default: 21)")
	umiLen := flag.Int("umi-len", 8, "UMI length (default: 8)")
	threads := flag.Int("t", defaultThreads(), "Worker goroutines")
	chunkSize := flag.Int("chunk-size", 200000, "Records per chunk")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Parallel BC+UMI extraction from R1 → R2 headers (umi_tools extract style)\n\n")
		fmt.Fprintf(os.Stderr, "Usage:\n")
		fmt.Fprintf(os.Stderr, "  umi_extract -1 r1.fq.gz -2 r2.fq.gz -o out.fq.gz [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		flag.PrintDefaults()
		fmt.Fprintf(os.Stderr, "\nExamples:\n")
		fmt.Fprintf(os.Stderr, "  # Default: bc 1-20bp, umi 21-28bp\n")
		fmt.Fprintf(os.Stderr, "  umi_extract -1 r1.fq.gz -2 r2.fq.gz -o out.fq.gz\n")
		fmt.Fprintf(os.Stderr, "  # Custom regions, 16 threads\n")
		fmt.Fprintf(os.Stderr, "  umi_extract -1 r1.fq -2 r2.fq -o out.fq \\\n")
		fmt.Fprintf(os.Stderr, "      --bc-start 1 --bc-len 16 --umi-start 17 --umi-len 10 -t 16\n")
	}
	flag.Parse()

	// Validate
	if *r1Path == "" || *r2Path == "" || *outPath == "" {
		fmt.Fprintln(os.Stderr, "[ERROR] -1, -2, -o are required")
		flag.Usage()
		os.Exit(1)
	}
	if *bcStart < 1 {
		fmt.Fprintln(os.Stderr, "[ERROR] --bc-start must be >= 1")
		os.Exit(1)
	}
	if *umiStart < *bcStart+*bcLen {
		fmt.Fprintf(os.Stderr, "[ERROR] --umi-start (%d) must >= bc-start + bc-len (%d)\n",
			*umiStart, *bcStart+*bcLen)
		os.Exit(1)
	}
	if *threads < 1 {
		*threads = 1
	}
	if *chunkSize < 1 {
		*chunkSize = 200000
	}
	for _, p := range []string{*r1Path, *r2Path} {
		if _, err := os.Stat(p); os.IsNotExist(err) {
			fmt.Fprintf(os.Stderr, "[ERROR] File not found: %s\n", p)
			os.Exit(1)
		}
	}

	fmt.Fprintf(os.Stderr, "[INFO] BC=%d-%d, UMI=%d-%d, threads=%d, chunk=%d\n",
		*bcStart, *bcStart+*bcLen-1,
		*umiStart, *umiStart+*umiLen-1,
		*threads, *chunkSize)

	// Open input files
	f1, err := openReader(*r1Path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "[ERROR] opening R1: %v\n", err)
		os.Exit(1)
	}
	defer f1.Close()

	f2, err := openReader(*r2Path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "[ERROR] opening R2: %v\n", err)
		os.Exit(1)
	}
	defer f2.Close()

	r1 := bufio.NewReaderSize(f1, 4*1024*1024)
	r2 := bufio.NewReaderSize(f2, 4*1024*1024)

	bw, cleanup, err := openWriter(*outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "[ERROR] opening output: %v\n", err)
		os.Exit(1)
	}
	defer cleanup()

	// ─── Pipeline: reader goroutine → worker pool → ordered writer ───

	chunkCh := make(chan Chunk, *threads*2)
	resultCh := make(chan Result, *threads*2)
	errCh := make(chan error, 1)

	// Reader goroutine
	go func() {
		chunkIdx := 0
		for {
			records, eof, rErr := readChunk(r1, r2, *chunkSize)
			if rErr != nil {
				errCh <- rErr
				close(chunkCh)
				return
			}
			if len(records) > 0 {
				chunkCh <- Chunk{Index: chunkIdx, Records: records}
				chunkIdx++
			}
			if eof {
				close(chunkCh)
				return
			}
		}
	}()

	// Worker goroutines
	var wg sync.WaitGroup
	for i := 0; i < *threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chunk := range chunkCh {
				res := processChunk(chunk, *bcStart, *bcLen, *umiStart, *umiLen)
				if res.Err != nil {
					errCh <- res.Err
					return
				}
				resultCh <- res
			}
		}()
	}

	// Close resultCh after all workers finish
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// Main: ordered writer
	pending := make(map[int]Result)
	nextIdx := 0
	totalWritten := 0
	tStart := time.Now()

	writeResult := func(res Result) {
		for _, rec := range res.Records {
			bw.WriteString(rec.Header)
			bw.WriteByte('\n')
			bw.WriteString(rec.Seq)
			bw.WriteByte('\n')
			bw.WriteString(rec.Plus)
			bw.WriteByte('\n')
			bw.WriteString(rec.Qual)
			bw.WriteByte('\n')
		}
		totalWritten += len(res.Records)
	}

mainLoop:
	for {
		select {
		case err := <-errCh:
			fmt.Fprintf(os.Stderr, "\n[ERROR] %v\n", err)
			os.Exit(1)

		case res, ok := <-resultCh:
			if !ok {
				break mainLoop
			}

			if res.Index == nextIdx {
				// Write this result, then flush any pending results now in order
				writeResult(res)
				nextIdx++
				for {
					if r, found := pending[nextIdx]; found {
						writeResult(r)
						delete(pending, nextIdx)
						nextIdx++
					} else {
						break
					}
				}
				elapsed := time.Since(tStart).Seconds()
				fmt.Fprintf(os.Stderr, "\r[MAIN] chunks_done=%d, records=%d, elapsed=%.1fs, %.0f records/s",
					nextIdx, totalWritten, elapsed, float64(totalWritten)/elapsed)
			} else {
				pending[res.Index] = res
			}
		}
	}

	elapsed := time.Since(tStart).Seconds()
	fmt.Fprintf(os.Stderr, "\n[INFO] Done! %d records written to %s, total time=%.1fs, %.0f records/s\n",
		totalWritten, *outPath, elapsed, float64(totalWritten)/elapsed)
}

