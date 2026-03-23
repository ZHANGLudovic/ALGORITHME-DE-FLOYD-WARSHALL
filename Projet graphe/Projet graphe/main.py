"""
Algorithme de Floyd-Warshall
SM601 - Théorie des graphes - Année 2025/2026
"""

import math
import sys
from pathlib import Path

INF = math.inf
BASE_DIR = Path(__file__).resolve().parent


# ──────────────────────────────────────────────
# 1. LECTURE DU GRAPHE DEPUIS UN FICHIER
# ──────────────────────────────────────────────

def lire_graphe(numero):
    """
    Lit le fichier graphe{numero}.txt et retourne :
      - n       : nombre de sommets
      - arcs    : liste de triplets (i, j, valeur)
    """
    nom_fichier = BASE_DIR / "graphes" / f"graphe{numero}.txt"
    with open(nom_fichier, "r") as f:
        lignes = [l.strip() for l in f if l.strip() != ""]

    n = int(lignes[0])
    nb_arcs = int(lignes[1])
    arcs = []
    for i in range(2, 2 + nb_arcs):
        parts = lignes[i].split()
        src, dst, val = int(parts[0]), int(parts[1]), int(parts[2])
        arcs.append((src, dst, val))

    return n, arcs


# ──────────────────────────────────────────────
# 2. CONSTRUCTION DE LA MATRICE D'ADJACENCE
# ──────────────────────────────────────────────

def construire_matrice(n, arcs):
    """
        Construit la matrice initiale des distances dist0 (n×n) :
            - dist0[i][i] = 0
            - dist0[i][j] = valeur de l'arc (i,j)  s'il existe
            - dist0[i][j] = +∞                    sinon
    """
    L = [[INF] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = 0
    for (src, dst, val) in arcs:
        L[src][dst] = val
    return L


# ──────────────────────────────────────────────
# 3. AFFICHAGE D'UNE MATRICE
# ──────────────────────────────────────────────

def afficher_matrice(M, n, titre=""):
    """Affiche une matrice de distances avec un rendu tabulaire aligné."""
    contenu = []
    for i in range(n):
        ligne = []
        for j in range(n):
            ligne.append("INF" if M[i][j] == INF else str(M[i][j]))
        contenu.append(ligne)
    afficher_tableau(contenu, n, titre=titre)


def afficher_matrice_P(P, n, titre=""):
    """Affiche la matrice suiv (sommet suivant) avec un rendu tabulaire aligné."""
    contenu = []
    for i in range(n):
        ligne = []
        for j in range(n):
            ligne.append("-" if P[i][j] is None else str(P[i][j]))
        contenu.append(ligne)
    afficher_tableau(contenu, n, titre=titre)


def afficher_tableau(contenu, n, titre=""):
    """Affiche un tableau n×n avec en-têtes de sommets (0..n-1) bien alignés."""
    if titre:
        print(f"\n  [  {titre}  ]")

    # Création des en-têtes de colonnes et de lignes explicitement
    entetes_colonnes = []
    for j in range(n):
        entetes_colonnes.append(str(j))

    entetes_lignes = []
    for i in range(n):
        entetes_lignes.append(str(i))

    # Détermination de la largeur des colonnes et de la colonne des indices
    largeur_col = 4
    for titre_col in entetes_colonnes:
        if len(titre_col) > largeur_col:
            largeur_col = len(titre_col)

    for ligne in contenu:
        for cellule in ligne:
            if len(cellule) > largeur_col:
                largeur_col = len(cellule)

    largeur_indice = 0
    for indice in entetes_lignes:
        if len(indice) > largeur_indice:
            largeur_indice = len(indice)

    # Ligne d'en-tête
    ligne_entete = "  " + " " * largeur_indice + " | "
    for index, col in enumerate(entetes_colonnes):
        if index > 0:
            ligne_entete += " | "
        ligne_entete += f"{col:>{largeur_col}}"
    print(ligne_entete)

    # Séparateur
    separator_lignes = "  " + "-" * largeur_indice + "-+-"
    for j in range(n):
        if j > 0:
            separator_lignes += "-+-"
        separator_lignes += "-" * largeur_col
    print(separator_lignes)

    # Corps du tableau
    for i in range(n):
        ligne_donnees = contenu[i]
        ligne = f"  {entetes_lignes[i]:>{largeur_indice}} | "
        for j in range(n):
            if j > 0:
                ligne += " | "
            ligne += f"{ligne_donnees[j]:>{largeur_col}}"
        print(ligne)

    print()


# ──────────────────────────────────────────────
# 4. ALGORITHME DE FLOYD-WARSHALL
# ──────────────────────────────────────────────

def copier_matrice(M):
    """Crée une copie profonde d'une matrice."""
    return [row[:] for row in M]


def floyd_warshall(L0, n):
    """
    Exécute l'algorithme de Floyd-Warshall.

    Retourne :
            - L_final : matrice finale des distances minimales
            - P_final : matrice finale suiv (sommet suivant)
            - historique_L : historique des matrices de distances
                                            (historique_L[0] = dist0, historique_L[1] = dist1, ...)
            - historique_P : historique des matrices suiv
                                            (historique_P[0] = suiv0, historique_P[1] = suiv1, ...)
    """
    L = copier_matrice(L0)

    # Initialisation de suiv0 : sommet suivant j si arc i->j, sinon None
    P = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and L[i][j] != INF:
                P[i][j] = j

    historique_L = [copier_matrice(L)]
    historique_P = [copier_matrice(P)]

    for k in range(n):
        for i in range(n):
            for j in range(n):
                if L[i][k] != INF and L[k][j] != INF:
                    nouveau = L[i][k] + L[k][j]
                    if nouveau < L[i][j]:
                        L[i][j] = nouveau
                        P[i][j] = P[i][k]
        historique_L.append(copier_matrice(L))
        historique_P.append(copier_matrice(P))

    return L, P, historique_L, historique_P


# ──────────────────────────────────────────────
# 5. DÉTECTION DES CIRCUITS ABSORBANTS
# ──────────────────────────────────────────────

def a_circuit_absorbant(L, n):
    """Un circuit absorbant existe si L[i][i] < 0 pour un certain i."""
    for i in range(n):
        if L[i][i] < 0:
            return True
    return False


# ──────────────────────────────────────────────
# 6. RECONSTITUTION ET AFFICHAGE D'UN CHEMIN
# ──────────────────────────────────────────────

def reconstruire_chemin(P, src, dst):
    """
    Reconstruit le chemin de src à dst via la matrice suiv.
    Retourne la liste des sommets du chemin, ou None si inexistant.
    """
    if P[src][dst] is None and src != dst:
        return None  # pas de chemin

    chemin = [src]
    courant = src
    # Garde-fou contre une matrice incohérente qui créerait une boucle.
    for _ in range(n := len(P)):
        if courant == dst:
            return chemin
        suivant = P[courant][dst]
        if suivant is None:
            return None  # chemin brisé
        courant = suivant
        chemin.append(courant)

    return None


def afficher_chemin(L, P, n, src, dst):
    """Affiche le chemin le plus court entre src et dst."""
    # Validation des sommets
    if not (0 <= src < n and 0 <= dst < n):
        print(f"    [!] Sommet invalide (0 a {n-1}).")
        return

    if src == dst:
        print(f"    [!] Depart et arrivee identiques ({src}).")
        return

    # Vérification de l'accessibilité
    if L[src][dst] == INF:
        print(f"    [!] Pas de chemin de {src} a {dst}.")
        return

    # Reconstruction et affichage du chemin
    chemin = reconstruire_chemin(P, src, dst)
    if chemin is None:
        print(f"    [!] Impossible de reconstruire le chemin.")
        return

    chemin_str = " -> ".join(str(s) for s in chemin)
    print(f"    [+] Chemin  : {chemin_str}")
    print(f"    [+] Valeur  : {L[src][dst]}")


# ──────────────────────────────────────────────
# 7. DESCRIPTIONS DES CAS
# ──────────────────────────────────────────────
DESCRIPTIONS_GRAPHES = {
    1: "Cas simple : 4 sommets, 5 arcs. Test basique de l'algorithme.",
    2: "Cas complexe : 5 sommets, 9 arcs. Arcs avec poids negatifs. Test general.",
    3: "Cas critique : 3 sommets, 4 arcs. Circuit absorbant detecte.",
    4: "Cas complet : 8 sommets, 16 arcs. Graphe de taille moyenne.",
    5: "4 sommets, 7 arcs. Graphe avec retour 2->0 et cycle 3<->2. Aucun circuit absorbant.",
    6: "4 sommets, 8 arcs. Boucle negative sur le sommet 1 (valeur -1). CIRCUIT ABSORBANT immediat.",
    7: "7 sommets, 12 arcs. Deux composantes non connectees (s0-s3 et s4-s6). Chemins infinis entre composantes.",
    8: "5 sommets, 7 arcs. Graphe acyclique (DAG), poids positifs. Plusieurs chemins alternatifs vers s4.",
    9: "5 sommets, 7 arcs. Arc negatif (2->3 = -2) et retour 4->1. Aucun circuit absorbant.",
    10: "5 sommets, 7 arcs. Cycle 2->3->4->2 de valeur totale -5 : CIRCUIT ABSORBANT.",
    11: "5 sommets, 8 arcs. Boucle positive sur s1 (valeur +2), arcs negatifs. Aucun circuit absorbant.",
    12: "5 sommets, 9 arcs. Arcs bidirectionnels 2<->3, arcs negatifs depuis s4. Aucun circuit absorbant.",
    13: "8 sommets, 12 arcs. Deux composantes non connectees (s0-s4 et s5-s7). Arcs negatifs. Cycle s5->s6->s5 = -1 : CIRCUIT ABSORBANT.",
    14: "10 sommets, 18 arcs. EXEMPLE APPLICATIF REEL : Reseau de transport urbain .",
}


def obtenir_description(numero):
    """Retourne la description du graphe ou un message par defaut."""
    return DESCRIPTIONS_GRAPHES.get(numero, "Graphe inconnu.")


# ──────────────────────────────────────────────
# 8. TRAITEMENT COMPLET D'UN GRAPHE
# ──────────────────────────────────────────────

def traiter_graphe(numero):
    print(f"\n{'='*64}")
    print(f"  GRAPHE N. {numero}")
    print(f"{'='*64}")
    
    # Description du cas
    description = obtenir_description(numero)
    print(f"\n  [INFO] {description}")

    # Lecture du graphe (on suppose pas d'erreur de format ou de fichier)
    n, arcs = lire_graphe(numero)
    print(f"  [+] Graphe charge : {n} sommet(s), {len(arcs)} arc(s)")

    # Construction de dist0 puis exécution de Floyd-Warshall
    print("  [*] Initialisation des matrices...\n")
    dist0 = construire_matrice(n, arcs)
    afficher_matrice(dist0, n, titre="dist0")
    
    print("  [*] Execution de l'algorithme de Floyd-Warshall...\n")
    L_final, P_final, hist_L, hist_P = floyd_warshall(dist0, n)

    # Affichage de suiv0 initialement
    afficher_matrice_P(hist_P[0], n, titre="suiv0")

    # Etapes successives
    for k in range(n):
        idx = k + 1
        afficher_matrice(hist_L[idx], n, titre=f"k={k} dist{idx}")
        afficher_matrice_P(hist_P[idx], n, titre=f"k={k} suiv{idx}")

    # Vérification des circuits absorbants
    if a_circuit_absorbant(L_final, n):
        print("  [!] CIRCUIT ABSORBANT DETECTE")
        print("  [!] La diagonale contient une valeur negative.")
        print("  [!] Les chemins minimaux ne peuvent pas etre determines.\n")
        return

    print("  [OK] Aucun circuit absorbant.\n")

    # Consultation interactive des chemins
    while True:
        rep = lire_chaine("  Afficher un chemin ? (o/n) : ")
        if rep not in ['o', 'oui']:
            break

        src = int(lire_chaine(f"    Sommet de depart  (0 a {n-1}) : "))
        dst = int(lire_chaine(f"    Sommet d'arrivee  (0 a {n-1}) : "))

        afficher_chemin(L_final, P_final, n, src, dst)
        print()


# ──────────────────────────────────────────────
# 9. BOUCLE PRINCIPALE
# ──────────────────────────────────────────────

def lire_chaine(prompt):
    """Lit une chaîne d'entrée de l'utilisateur.

    Si l'entrée n'est pas interactive (mode pipe/fichier), renvoie une chaîne vide.
    """
    if not sys.stdin.isatty():
        return ""
    return input(prompt).strip().lower()


def main():
    
    print("|                                                              |")
    print("|        Algorithme de FLOYD-WARSHALL                          |")
    print("|        Recherche des chemins les plus courts                 |")
    print("|                                                              |")
    

    while True:
        print("-" * 64)
        rep = lire_chaine("  Traiter un graphe ? (o/n) : ")
        if rep not in ['o', 'oui']:
            print("\n  A bientot!\n")
            break

        num = int(lire_chaine("  Numero du graphe (1-14) : "))
        # On suppose que l'utilisateur saisit un entier valide entre 1 et 13

        traiter_graphe(num)


if __name__ == "__main__":
    main()