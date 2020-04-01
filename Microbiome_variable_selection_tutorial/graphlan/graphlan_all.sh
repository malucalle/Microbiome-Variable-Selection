python ./graphlan_functions/graphlan_annotate.py --annot annot_0.txt taxa.txt taxa_1.xml
python ./graphlan_functions/graphlan_annotate.py --annot annot_all.txt taxa_1.xml taxa_2.xml
python ./graphlan_functions/graphlan.py taxa_2.xml taxa.png --dpi 300 --size 3.5
python ./graphlan_functions/graphlan.py taxa_2.xml taxa.svg --dpi 300 --size 3.5
