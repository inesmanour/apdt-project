# SUJET : ASSIGNATION ET DETECTION DES PARTIES TRANSMEMBRANAIRES D'UNE PROTÉINE

**Auteur :** Inès MANOUR  
**Niveau :** Master 2 Bioinformatique  
**Année :** 2025-2026

---

## 🧬 Contexte

Les protéines transmembranaires jouent un rôle clé dans la signalisation cellulaire et les échanges moléculaires.  
L’objectif de ce projet est d’**implémenter une méthode automatique** pour :

- Identifier les régions transmembranaires à partir de structures PDB.  
- Associer hydrophobicité et accessibilité au solvant (RSA, via DSSP).  
- Estimer la position et l’orientation de la membrane.  
- Visualiser la membrane et les dalles associées en **PyMOL** ou **Matplotlib**.  

La méthodologie suit la publication de **Tusnády–Dosztányi–Simon (2004, Bioinformatics)** et utilise la **scale de Kyte–Doolittle** pondérée par RSA.

---

## 📂 Organisation du dépôt 

```
adpt-project/
├── data/                     
│   ├── 1K24.pdb
│   └── 1K24.dssp
│
├── results/             
│   ├── 2025-09-10/      
│   │   ├── 1K24_result.json
│   │   ├── 1K24_annot.csv
│   │   └── 1K24_view.png
│   
│
├── src/                 # Code source Python
│   ├── main.py          # Script principal (CLI)
│   ├── protein.py       # Parsing PDB/DSSP, extraction Cα
│   ├── geometrie.py     # Projection, slabs, optimisation
│   ├── hydrophobicite.py # Scores hydrophobes (KD × RSA)
│   ├── visualisation.py # Matplotlib + PyMOL .pml
│   └── init.py
│   └─ envi.yml 
│
└── README.md            
```

---

## ⚙️ Installation

### 1. Cloner le dépôt
```
git clone https://github.com/ines-manour/adpt-project.git
cd adpt-project
```


### 2. Créer un environnement Python

#### Option — avec Pixi 

```
pixi init 
pixi install
pixi add python>=3.10 numpy matplotlib biopython dssp
pixi add pymol-open-source
pixi shell
```


### 2. Lancer l’analyse

python src/main.py \
  --pdb data/1K24.pdb \
  --dssp data/1K24.dssp \
  --n_dirs 30 \
  --head 5 \
  --plot \
  --widen-step 1 \
  --max-half-extra 30 \
  --pml view_membrane.pml


### 3. Visualiser les résultats

#### Exemple avec 1k24pdb :  

```
python main.py --pdb data/1k24.pdb --n_dirs 30 --head 5 --plot   --widen-step 1 --max-half-extra 30 --pml view_membrane.pml
```

```
  •	--pdb → la protéine à analyser
	•	--n_dirs → nombre de directions testées
	•	--head → combien de résidus afficher avec ASA/RSA
	•	--plot → afficher ou non la figure matplotlib
	•	--widen-step → pas d’élargissement de la dalle
	•	--max-half-extra → jusqu’où on peut élargir
	•	--pml → exporter un script pour PyMOL
```

##### Résultats : 

##### • En console :

    ```
    === 1k24.pdb ===
    - Cα lus : 249
    - ASA (DSSP) : 249/249
    - RSA moyenne : 0.288

    --- Premiers résidus (ASA/RSA) ---
    A (' ', 5, ' ') GLN | ASA=235.0 | RSA=1.000
    A (' ', 6, ' ') THR | ASA=73.0 | RSA=0.424
    A (' ', 7, ' ') ALA | ASA=65.0 | RSA=0.504
    A (' ', 8, ' ') ASN | ASA=82.0 | RSA=0.421
    A (' ', 9, ' ') GLU | ASA=52.0 | RSA=0.233
    - Centre de masse (Cα) : 30.403, -2.402, 9.879

    --- Scores par direction (Q) ---
    1. Q=0.731 | n=[0.256 0.    0.967] | centre=0.26 Å | demi-ép=12.50 Å | d_in=-12.24 Å | d_out=12.76 Å
    2. Q=0.694 | n=[ 0.227  0.119 -0.967] | centre=0.23 Å | demi-ép=13.50 Å | d_in=-13.27 Å | d_out=13.73 Å
    3. Q=0.563 | n=[-0.86   0.355  0.367] | centre=-0.86 Å | demi-ép=28.50 Å | d_in=-29.36 Å | d_out=27.64 Å
    4. Q=0.503 | n=[ 0.048 -0.551  0.833] | centre=0.05 Å | demi-ép=9.50 Å | d_in=-9.45 Å | d_out=9.55 Å
    (...)

    → Q_max = 0.731
    - Normale estimée              : [0.25603819 0.         0.96666667]
    - Centre optimal (le long de n): 0.26 Å
    - Score max (KD×RSA moy.)      : 0.731
    - Demi-épaisseur optimale      : 12.50 Å (épaisseur ≈ 25.00 Å)
    - Positions d_in / d_out       : -12.24 Å / 12.76 Å le long de n
    - Centres des dalles (3D)      : IN  = [27.268, -2.402, -1.957] | OUT = [33.669, -2.402, 22.210]
    [PyMOL] Script écrit → view_membrane.pml  | Dans PyMOL : @view_membrane.pml
    - Script PyMOL écrit : view_membrane.pml
    ```


##### •	En PyMOL :
    
**Figure PyMOL :**
![Membrane dans PyMOL](data/my_membrane.png)


•	En Matplotlib :
Si --plot est activé, un rendu 3D interactif s’ouvre.

**Figure Matplotlib :**
![Membrane dans Matplotlib](data/Figure_1plot.png)




## 🎓 Notes pédagogiques

	•	Le code est écrit selon PEP 8 et documenté selon PEP 257.
	•	Usage de @dataclass pour simplifier les classes de données (ResidueCA).
	•	Modularité : chaque fichier Python est spécialisé (proteine, géométrie, hydrophobicité, visualisation).
	•	Paramètres ajustables : seuil RSA, épaisseur, nombre de directions.
	


## 📚 Références
	•	Tusnády, G.E., Dosztányi, Z., Simon, I. (2004). Transmembrane proteins in the Protein Data Bank: identification and classification. Bioinformatics, 20(17), 2964–2972.
	•	Kyte, J., & Doolittle, R.F. (1982). A simple method for displaying the hydropathic character of a protein. J Mol Biol, 157(1), 105–132.
	•	DSSP: Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure. Biopolymers, 22(12), 2577–2637.