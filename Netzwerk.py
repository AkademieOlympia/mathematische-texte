import networkx as nx

# Initialisierung des Graphen für die Netzwerkanalyse
G = nx.Graph()

# Definition der zentralen Knoten (Hubs) und Brücken (Bridges)
# Basierend auf den Epstein-Akten und den Homburg-Leaks
connections = [
    ("Jeffrey Epstein", "Lesley Groff"),
    ("Jeffrey Epstein", "Melanie Walker"),
    ("Lesley Groff", "Philippa Sigl-Glöckner"), # Leak: Matrixfehler
    ("Melanie Walker", "Philippa Sigl-Glöckner"), # Leak: Starfriseur-Kontakt
    ("Jeffrey Epstein", "Nicole Junkermann"),
    ("Nicole Junkermann", "Peter Thiel"),
    ("Nicole Junkermann", "Ehud Barak"),
    ("Ehud Barak", "Carbyne/Reporty"),
    ("Jeffrey Epstein", "Carbyne/Reporty"),
    ("Philippa Sigl-Glöckner", "BMF/Olaf Scholz"), # Politische Platzierung
    ("Nicole Junkermann", "Infront/Adidas"),
]

G.add_edges_from(connections)

def analyze_network(graph):
    # Berechnung der Zentralitätswerte
    # Betweenness Centrality zeigt die Rolle als Vermittler/Gatekeeper
    betweenness = nx.betweenness_centrality(graph)
    
    # Clustering Coefficient zeigt, wie stark die Umgebung vernetzt ist
    clustering = nx.clustering(graph)
    
    print(f"{'Akteur':<25} | {'Betweenness':<12} | {'Clustering'}")
    print("-" * 55)
    for node in sorted(betweenness, key=betweenness.get, reverse=True):
        print(f"{node:<25} | {betweenness[node]:.4f} | {clustering[node]:.4f}")

# Analyse ausführen
analyze_network(G)