#!/usr/bin/env python3
"""
Obsidian Vault Gap Analyzer & Link Suggester
--------------------------------------------
Analysiert einen Obsidian-Vault auf unentdeckte Wissenslücken, unvollständige Verknüpfungen 
(Missing Links) und isolierte Themeninseln mittels TF-IDF / Cosine Similarity.
"""

import os
import re
import json
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Set
import numpy as np

try:
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.metrics.pairwise import cosine_similarity
except ImportError:
    TfidfVectorizer = None
    cosine_similarity = None

class ObsidianVaultAnalyzer:
    def __init__(self, vault_path: str):
        self.vault_path = Path(vault_path).expanduser()
        self.notes: Dict[str, Dict] = {}
        self.titles: List[str] = []
        self.paths: List[Path] = []
        self.contents: List[str] = []
        
    def scan_vault(self):
        """Liest alle Markdown-Dateien im Vault ein und extrahiert Links und Tags."""
        if not self.vault_path.exists():
            raise FileNotFoundError(f"Vault-Pfad {self.vault_path} existiert nicht.")
            
        for path in self.vault_path.rglob("*.md"):
            try:
                content = path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
                
            title = path.stem
            
            # Extrahiere Wiki-Links [[link]] oder [[link|alias]]
            wiki_links = set(re.findall(r"\[\[([^\]\|]+)(?:\|[^\]]+)?\]\]", content))
            # Extrahiere Tags #tag
            tags = set(re.findall(r"#([a-zA-Z0-9_\-/]+)", content))
            
            # Säubere den Text für NLP (Entferne Codeblöcke und Links-Syntax)
            clean_text = re.sub(r"```.*?```", "", content, flags=re.DOTALL)
            clean_text = re.sub(r"\[\[.*?\]\]", "", clean_text)
            clean_text = re.sub(r"#\S+", "", clean_text)
            
            self.notes[title] = {
                "title": title,
                "path": path,
                "content": content,
                "clean_text": clean_text,
                "links": wiki_links,
                "tags": tags
            }
            
            self.titles.append(title)
            self.paths.append(path)
            self.contents.append(clean_text)
            
    def analyze_gaps(self, similarity_threshold: float = 0.25, top_k: int = 5) -> Dict:
        """Sucht nach semantisch ähnlichen Zetteln, die bisher noch NICHT miteinander verlinkt sind."""
        if not self.titles:
            return {"error": "Keine Notizen im Vault gefunden."}
            
        if TfidfVectorizer is None:
            return {"error": "scikit-learn ist nicht installiert."}
            
        vectorizer = TfidfVectorizer(stop_words='english', min_df=1)
        try:
            tfidf_matrix = vectorizer.fit_transform(self.contents)
        except ValueError:
            # Falls Texte zu kurz oder stumme Wörter
            return {"error": "Textmaterial im Vault unzureichend für TF-IDF."}
            
        sim_matrix = cosine_similarity(tfidf_matrix)
        
        missing_links = []
        orphan_notes = []
        cluster_suggestions = []
        
        n = len(self.titles)
        for i in range(n):
            title_a = self.titles[i]
            note_a = self.notes[title_a]
            existing_links = note_a["links"]
            
            # Prüfe auf verwaiste Notizen (Orphans)
            incoming_links = sum(1 for t, info in self.notes.items() if title_a in info["links"])
            outgoing_links = len(existing_links)
            
            if outgoing_links == 0 and incoming_links == 0:
                orphan_notes.append(title_a)
                
            for j in range(i + 1, n):
                title_b = self.titles[j]
                note_b = self.notes[title_b]
                
                sim_score = sim_matrix[i, j]
                
                # Ist bereits eine Verlinkung vorhanden (A -> B oder B -> A)?
                is_linked = (title_b in note_a["links"]) or (title_a in note_b["links"])
                
                if sim_score >= similarity_threshold and not is_linked:
                    missing_links.append({
                        "note_a": title_a,
                        "note_b": title_b,
                        "similarity": round(float(sim_score), 4),
                        "suggested_link": f"[[{title_b}]] in [[{title_a}]] einfügen"
                    })
                    
        # Sortiere fehlende Verknüpfungen nach Ähnlichkeit
        missing_links.sort(key=lambda x: x["similarity"], reverse=True)
        
        return {
            "total_notes": n,
            "orphan_notes": orphan_notes,
            "missing_links_count": len(missing_links),
            "top_missing_links": missing_links[:top_k * 3],
        }

    def generate_report_markdown(self, analysis_result: Dict) -> str:
        """Erzeugt einen detaillierten Markdown-Bericht zur Integration in Obsidian."""
        md = []
        md.append("# 🕸️ Obsidian Vault Knowledge Gap Report")
        md.append(f"**Analysierte Notizen:** {analysis_result.get('total_notes', 0)}\n")
        
        md.append("## 🔍 Empfohlene Querverbindungen (Missing Links)")
        md.append("Diese Notizen weisen eine hohe inhaltliche/semantische Ähnlichkeit auf, sind aber aktuell **nicht miteinander verlinkt**:\n")
        
        top_links = analysis_result.get("top_missing_links", [])
        if top_links:
            md.append("| Ähnlichkeit | Notiz A | Notiz B | Empfohlene Aktion |")
            md.append("|-------------|---------|---------|-------------------|")
            for item in top_links:
                md.append(f"| {item['similarity']*100:.1f}% | [[{item['note_a']}]] | [[{item['note_b']}]] | `[[{item['note_b']}]]` zu [[{item['note_a']}]] hinzufügen |")
        else:
            md.append("*(Keine unentdeckten Wissenslücken gefunden – dein Vault ist hervorragend vernetzt!)*\n")
            
        md.append("\n## 🏝️ Verwaiste Notizen (Orphan Notes)")
        md.append("Notizen ohne jegliche eingehende oder ausgehende Verlinkungen:\n")
        
        orphans = analysis_result.get("orphan_notes", [])
        if orphans:
            for orphan in orphans:
                md.append(f"- [[{orphan}]]")
        else:
            md.append("*(Keine verwaisten Notizen gefunden)*\n")
            
        return "\n".join(md)

if __name__ == "__main__":
    import sys
    
    vault = sys.argv[1] if len(sys.argv) > 1 else "./test_vault"
    print(f"Starte Wissenslücken-Analyse für Vault: {vault}")
    
    analyzer = ObsidianVaultAnalyzer(vault)
    try:
        analyzer.scan_vault()
        results = analyzer.analyze_gaps()
        report = analyzer.generate_report_markdown(results)
        
        report_file = Path(vault) / "00_Knowledge_Gap_Report.md"
        report_file.parent.mkdir(parents=True, exist_ok=True)
        report_file.write_text(report, encoding="utf-8")
        print(f"Bericht erfolgreich geschrieben: {report_file}")
    except Exception as e:
        print(f"Fehler bei der Analyse: {e}")
