Welcome to Multi-PredGO Application:
Multi-PredGO is a multi-modal protein function prediction model that uses the protein sequence, protein structure, and protein-protein interaction network-based information to predict GO-based protein function. As the protein function classes are dependent on each other, we have used a neuro-symbolic hierarchical classification model, which resembles the structure of Gene Ontology (GO), for effectively predicting the dependent protein functions.

Installation:
command 1 : conda create -n multi-predGO python=2.7
command 2 : source activate multi-predGO
command 3 : conda install -c bioconda diamond
command 4 : pip install -r requirements.txt
command 5 : pip install Flask
command 6 : pip install -U flask-cors



Run the code:
1) To run the code as application: 
Run app.py using command python app.py
app.py will run ur application  in url "http://0.0.0.0:5000/", you can change the host and the part address in app.py
you can get the predictions uisng  http://0.0.0.0:5000/accession?accession_no={accession_no}
example:  http://0.0.0.0:5000/accession?accession_no=P31946, here P31946 is the accession no corresponding which we want a prediction


