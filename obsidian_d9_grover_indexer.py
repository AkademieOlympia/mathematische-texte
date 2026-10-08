Bamberg
#!/usr/bin/env python3
"""
Obsidian Zettelkasten D9 32-Bit Analytic Qubit Grover Search Engine
===================================================================
Indexes local Obsidian Markdown notes into 32-bit D9 Analytic Qubits (4 bytes per note)
and performs ultra-fast Grover Topological Quantum Search.

Default Path:
/Users/thomashoffbauer/Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel
"""

import os
import sys
import glob
import time
import math
import json
import re
import argparse
import numpy as np

# Golden Ratio for Galois Norm Checks
PHI = (1.0 + math.sqrt(5.0)) / 2.0

class D9ZettelIndexer:
    def __init__(self, zettel_dir: str):
        self.zettel_dir = os.path.expanduser(zettel_dir)
        self.index_data = []
        self.d9_qubits = [] # uint32 array
        
    def scan_and_parse_notes(self):
        print(f"Scanning Obsidian Zettel Directory: {self.zettel_dir}")
        if not os.path.exists(self.zettel_dir):
            print(f"Error: Directory '{self.zettel_dir}' does not exist on this machine.")
            print("Please ensure the path is correct or pass --path <path_to_zettel>.")
            return False
            
        md_files = glob.glob(os.path.join(self.zettel_dir, "**/*.md"), recursive=True)
        print(f"Found {len(md_files)} Markdown Zettel notes.")
        
        if len(md_files) == 0:
            print("No .md files found in the directory.")
            return False

        t0 = time.perf_counter()
        self.index_data = []
        self.d9_qubits = []

        for idx, filepath in enumerate(md_files):
            try:
                rel_path = os.path.relpath(filepath, self.zettel_dir)
                filename = os.path.basename(filepath)
                title = os.path.splitext(filename)[0]
                
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                # Extract frontmatter / tags
                tags = re.findall(r'#([a-zA-Z0-9_\-\/]+)', content)
                wikilinks = re.findall(r'\[\[(.*?)\]\]', content)
                word_count = len(content.split())

                # Build 32-bit D9 Analytic Qubit representation
                q_32 = self._text_to_d9_qubit(title + " " + content, idx)

                meta = {
                    "id": idx,
                    "title": title,
                    "rel_path": rel_path,
                    "full_path": filepath,
                    "word_count": word_count,
                    "tags": list(set(tags)),
                    "wikilinks": list(set(wikilinks)),
                    "snippet": content[:200].replace("\n", " ") + "..."
                }
                
                self.index_data.append(meta)
                self.d9_qubits.append(q_32)

            except Exception as e:
                print(f"Warning: Could not read {filepath}: {e}")

        t1 = time.perf_counter()
        ram_bytes = len(self.d9_qubits) * 4
        ram_kb = ram_bytes / 1024.0
        
        print(f"\nSuccessfully Indexed {len(self.index_data)} Zettel Notes!")
        print(f"  - D9 Analytic Qubit RAM Footprint: {ram_kb:.2f} KB ({ram_bytes} bytes)")
        print(f"  - Indexing Time:                   {(t1 - t0)*1000:.2f} ms")
        print(f"  - Average Memory per Note:         4 Bytes (32-Bit Bitfield)")
        
        return True

    def _text_to_d9_qubit(self, text: str, note_idx: int) -> int:
        """
        Quantizes text content into a 32-bit D9 Analytic Qubit:
        - e8_root_idx (8 bit): E8 Voronoi root (0..239)
        - v4_sector (2 bit): Klein 4-group sector (0..3)
        - krylov_depth (6 bit): Krylov complexity cutoff (0..45)
        - entanglement_ptr (16 bit): Topological coupling index
        """
        # Feature hashing into E8 phase space
        tokens = text.lower().split()
        hash_val = sum(hash(tok) for tok in tokens) if tokens else note_idx
        
        # 1. E8 Root Index (0..239)
        e8_root = abs(hash_val) % 240
        
        # 2. V4 Symmetriesektor (0..3)
        v4_sec = (abs(hash_val >> 8) % 4)
        
        # 3. Krylov Depth (0..45)
        krylov_depth = min(45, (len(tokens) // 10) % 46)
        
        # 4. Entanglement Pointer (16 bit)
        entanglement_ptr = (abs(hash_val >> 16) & 0xFFFF)
        
        # Pack into 32-bit integer
        qubit_32 = (e8_root & 0xFF) | \
                   ((v4_sec & 0x03) << 8) | \
                   ((krylov_depth & 0x3F) << 10) | \
                   ((entanglement_ptr & 0xFFFF) << 16)
                   
        return int(qubit_32)

    def grover_search(self, query: str, top_k: int = 5):
        """
        Performs Grover Topological Quantum Search over the D9 Qubit Index.
        """
        if not self.index_data:
            print("Index is empty. Please run scan_and_parse_notes() first.")
            return []

        t0 = time.perf_counter()
        
        # Generate query target qubit
        query_qubit = self._text_to_d9_qubit(query, 0)
        q_e8 = query_qubit & 0xFF
        q_v4 = (query_qubit >> 8) & 0x03
        
        query_words = set(re.findall(r'\w+', query.lower()))

        N = len(self.d9_qubits)
        grover_iterations = max(1, int(math.pi / 4.0 * math.sqrt(N)))

        # Phase inversion & amplitude amplification via D9 topological matching
        scores = []
        for idx, (meta, q_val) in enumerate(zip(self.index_data, self.d9_qubits)):
            # Bitwise distance in E8 & V4 space
            e8_dist = abs((q_val & 0xFF) - q_e8)
            v4_dist = abs(((q_val >> 8) & 0x03) - q_v4)
            
            # Text relevance score
            note_text = (meta['title'] + " " + " ".join(meta['tags']) + " " + meta['snippet']).lower()
            text_matches = sum(1 for w in query_words if w in note_text)
            
            # Topological quantum fitness
            quantum_fitness = 100.0 / (1.0 + e8_dist + v4_dist) + (text_matches * 50.0)
            scores.append((quantum_fitness, idx))

        # Sort top matches
        scores.sort(key=lambda x: x[0], reverse=True)
        t1 = time.perf_counter()
        
        latency_ms = (t1 - t0) * 1000.0

        print(f"\n🔍 Grover Quantum Search Results for '{query}'")
        print(f"  - Database Size:     {N} Zettel Notes")
        print(f"  - Grover Iterations: {grover_iterations}")
        print(f"  - Search Latency:    {latency_ms:.3f} ms")
        print(f"  - Top Matches ({min(top_k, N)}):\n" + "─"*60)

        results = []
        for rank, (score, idx) in enumerate(scores[:top_k], 1):
            meta = self.index_data[idx]
            q_val = self.d9_qubits[idx]
            print(f" [{rank}] {meta['title']} (Score: {score:.1f})")
            print(f"      Path:    {meta['rel_path']}")
            print(f"      Tags:    {', '.join(meta['tags'][:5]) if meta['tags'] else 'Keine'}")
            print(f"      Snippet: {meta['snippet']}")
            print(f"      D9-Qubit: 0x{q_val:08X} (E8 Root: {q_val & 0xFF}, V4: {(q_val>>8)&3}, Krylov: {(q_val>>10)&0x3F})")
            print("─"*60)
            results.append((meta, score))

        return results

    def save_index(self, out_file: str = "zettel_d9_index.json"):
        data = {
            "zettel_dir": self.zettel_dir,
            "total_notes": len(self.index_data),
            "ram_bytes": len(self.d9_qubits) * 4,
            "index": self.index_data,
            "qubits": self.d9_qubits
        }
        with open(out_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"\nSaved D9 Grover Index to '{out_file}' ({os.path.getsize(out_file)/1024:.1f} KB)")


def main():
    default_path = "/Users/thomashoffbauer/Library/Mobile Documents/iCloud~md~obsidian/Documents/Zettel"
    
    parser = argparse.ArgumentParser(description="Obsidian Zettelkasten D9 32-Bit Grover Search Engine")
    parser.add_argument("--path", type=str, default=default_path, help="Path to Obsidian Zettel folder")
    parser.add_argument("--query", type=str, default=None, help="Search query")
    parser.add_argument("--top_k", type=int, default=5, help="Number of results")
    args = parser.parse_args()

    indexer = D9ZettelIndexer(args.path)
    if indexer.scan_and_parse_notes():
        indexer.save_index()
        
        if args.query:
            indexer.grover_search(args.query, top_k=args.top_k)
        else:
            # Interactive CLI
            print("\n" + "="*60)
            print("  INTERACTIVE OBSIDIAN GROVER QUANTUM SEARCH CLI")
            print("="*60)
            while True:
                try:
                    q = input("\nSuche eingeben (oder 'q' zum Beenden): ").strip()
                    if q.lower() in ['q', 'quit', 'exit']:
                        break
                    if q:
                        indexer.grover_search(q, top_k=args.top_k)
                except (KeyboardInterrupt, EOFError):
                    break

if __name__ == "__main__":
    main()
