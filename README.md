# projet-VDB

Evaluation VDB : 15 décembre 


_Dépot portable sur git : Notebook --> interacftif, figure a l'intérieur

_Poster : Clarté (+ ou - 4 figures pertinante, retravaillé, combiné) + originalité +  discours (corrélation, question scientifique simple, représentation graphique, exploration de donnée formuler des hypothèse scienditifque) + 

_ Optionnel (bonnus) : production d'une petite app shiny pour interactivité via seuil, profondeur,...


15 décembre : 11H30/15H --> présentation poster (10 minutes de présentation / 5 minutes de question)


DEMARCHE :
Etape 1 : Construire des visualisations simples
_On occulte l’information hiérarchique et on fixe un niveau taxonomique (e.g. Genre)
_Stackedbarplots
_heatmap


Etape 2 : Construire des visualisations plus complexes
_On exploite l’information hiérarchique
_Treemap
_Sunburst diagram (circular treemap)


Etape 3 : Interactivité
_Permettre à l’utilisateur de choisir l’abondance minimale(slider) pour représenter un OTU et afficher le stackedbarplotet la heatmap
_Permettre à l’utilisateur de choisir le niveau taxonomique(select box) pour afficher le stackedbarplotet la heatmap
_Permettre à l’utilisateur de rechercher un taxon (e.g. Escherichia) et afficher son abondance dans les différents échantillons sous forme de barplot
_Afficher des informations lors du survol de la souris•Générer des tableaux de données téléchargeable


Etape 4 : Extraire de l’information
_Chercher une corrélation entre les variables
_Type d’échantillon (sol, toilettes, mains)
_Profil taxonomique global
_Abondance d’un taxa
_Analyses multivariées (e.g. ACP)
_Etc.



METHODE :
_Jupyter Notebook : fichier ipynb
_Outils/Package python utilisable : Jupyter architecture & interface & widget API, Pandas, Scipy, Seaborn, matplotlib (figure and axis, colors, callbacks)

