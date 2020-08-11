# MultiPredGO

## This is the readme file that contains the guidelines and information about the compilation of the code of the following paper

**Paper Name:-** MultiPredGO: Deep Multi-Modal Protein Function Prediction by Integrating Protein Structure, Sequence, and Interaction Information

- **Authors:** Swagarika Jaharlal Giri<sup>1</sup>, Pratik Dutta<sup>1</sup>, Parth Halani<sup>2</sup> and Sriparna Saha<sup>1</sup>
- **Affiliation:** <sup>1</sup>Department of Computer Science and Engineering, IIT Patna, India, <sup>2</sup>Depatment of Computer Science and Engineering, IIIT Guwahati, India
- **Web-server url 1:** http://multipred.co.in
- **Web-server url 2:** http://13.233.6.192:3000/
- **Web-server url 3:** http://www1.iitp.ac.in/~pratik.pcs16/multipred.html
- **Web-server Backend API:** http://15.206.210.184/accession?accession_no={Uniprot's Accession No}&ontology={ontology}
- **Web-server Backend API example:** http://15.206.210.184/accession?accession_no=P31946&ontology=cc

- **Status:** Under Minnor Revison in [IEEE Journal of Biomedical and Health Informatics(IEEE JBHI)](https://jbhi.embs.org/)

## Welcome to MultiPredGO Application:

MultiPredGO is a multi-modal protein function prediction model that uses the protein sequence, protein structure, and protein-protein interaction network-based information to predict GO-based protein function. As the protein function classes are dependent on each other, we have used a neuro-symbolic hierarchical classification model, which resembles the structure of Gene Ontology (GO), for effectively predicting the dependent protein functions.

## Installation:

```bash
conda create -n multi-pred python=2.7
source activate multi-pred
conda install -c bioconda diamond
pip install -r requirements.txt
pip install Flask
pip install -U flask-cors
```

## File Organization

```
|---- data
|    |--- models
|         | ---  model_bp.h5      <----- the  best saved model for BP (Multi-predGO)
|         | ---  model_cc.h5      <----- the  best saved model for CC(Multi-predGO)
|         | ---  model_mf.h5      <----- the  best saved model for MF(Multi-predGO)
|         | ---  mode_seq_bp.h5   <----- the  best saved model for Sequential BP
|         | ---  mode_seq_cc.h5   <----- the  best saved model for Sequential CC
|         | ---  mode_seq_mf.h5   <----- the  best saved model for Sequential MF
|
|
|    |---        AccessionNumber_Structure_StatusFileWithAccessionIndex.pkl <-  the file
|                                                         that has the status
|                                                         of accession no. status
|
|    |---        AccessionNumberStatusFileWithAccessionIndex.pkl <- "the file that has the
|                                                           status of accession no. if the
|  accession number's sequence , structure and PPIN feature is present in our database or not.
|
|    |---       bp.pkl   <--- The sequence of GO term/ Function sequence (BP)
|    |---       cc.pkl   <--- The sequence of GO term/ Function sequence (CC)
|    |---       mf.pkl   <--- The sequence of GO term/ Function sequence (MF)
|
|    |---       combined-multimodal-bp.pkl   <--- this dataset contains the fasta sequence,
|                   PPIN Knowledge Graph based embedding, and Structural feature information for proteins
|                   participating in biological process.
|
|    |---        combined-multimodal-mf.pkl   <--- this dataset contains the fasta sequence,
|                   PPIN Knowledge Graph based embedding, and Structural feature information for proteins
|                   participating in molecular function.
|
|    |---       combined-multimodal-bp.pkl   <--- this dataset contains the fasta sequence,
|                  PPIN Knowledge Graph based embedding, and Structural feature information for proteins
|                  participating in biological process.
|    |---         go.obo     <--- gene ontology
|    |---         go.txt     <--- gene ontology text file
|    |---         interpro.xml <--- file used to make PPI interection network
|
|----  __init__.py  <--- initiate file
|
|----  app.py      <---- the main file that we need to run using command  (python app.py)
|                         port configuration and URL.
|
|----  plots-folder  <----  plot the roc curve for the following predictions
|
|----  DataGenerationAPI.py    <---- It is the file that is called by app.py handles creation of json object after concatenating all the results
|----  ForTestMultiPredModelAPI.py  <---- It is the file where model is loaded and then results are calculated
|----  giri_inga_formatted   <---- prediction of functions(GO terms) after applying INGA model
|----  inga-predictions-bp.csv  <---- inga-prediction file for Biological Process
|----  inga-predictions-cc.csv  <---- inga-prediction file for Cellular Component
|----  inga-predictions-mf.csv  <---- inga-prediction file for Molecular Function
|----  pred_bp.pkl  <---- Multi-PredGO prediction(BP)
|----  pred_mf.pkl  <---- Multi-PredGO prediction(MF)
|----  pred_cc.pkl  <---- Multi-PredGO prediction(CC)
|----  utils.py    <---- has all utils functions of creating ontology etc
|----  List_Accession_No.txt   <----- contains all the accession no for which we were able to get structure, sequence and PPIN information and
|                                     the  features are stored.
|
|
```

## Run the code:

```
## Run as a API in linux/ubuntu/macOS/Windows

****** Approach 1 : Run as a API in linux/ubuntu/macOS/Windows ***************************************
______________________________________________________________________________________________________
Step 1: Installation
______________________________________________________________________________________________________
conda create -n multi-predGO python=2.7
source activate multi-predGO
conda install -c bioconda diamond
pip install -r requirements.txt
pip install Flask
pip install -U flask-cors
_______________________________________________________________________________________________________
Step 2 : Run the application

python app.py   <--------- Runs the application on  "http://0.0.0.0:5000/"

_______________________________________________________________________________________________________
Step 3 : Get Results for a single protein (You can choose any accession no from List_Accession_No.txt)

## using any Web browser
Type the  url
http://0.0.0.0:5000/accession?accession_no={accession_no}&ontology={ontology}
ex: http://0.0.0.0:5000/accession?accession_no=P31946&ontology=cc


______________________________________________________________________________________________________

## Get Results through a script

************************ Approach 2 : Run as script in linux/ubuntu/macOS/Windows ********************

______________________________________________________________________________________________________
Step 1: Installation
______________________________________________________________________________________________________
conda create -n multi-predGO python=2.7
source activate multi-predGO
conda install -c bioconda diamond
pip install -r requirements.txt
pip install Flask
pip install -U flask-cors
______________________________________________________________________________________________________

step 2:
 > python DataGenerationAPI.py       <------ run this DataGenerationAPI script
 > Please enter the accession no P31946. <----- Can use any accession no corresponding to protein
 > Please enter the ontology bp.  <------ Desired Ontology for which results are to be received


_____________________________________________________________________________________________________

 ## Train Model and Test the model
 ************************ Approach 3 : Run as script in linux/ubuntu/macOS/Windows(Train/Test) *********
 ______________________________________________________________________________________________________
Step 1: Installation
______________________________________________________________________________________________________
conda create -n multi-predGO python=2.7
source activate multi-predGO
conda install -c bioconda diamond
pip install -r requirements.txt

______________________________________________________________________________________________________

step 2: Run the script

python MultiPredGO.py.  <-----  script that takes protein sequence, structure PPIN based knowledge
                                  graph embedding to predict the protein function(GO term).
                                  uses data/multimodaltrain-{ontology}.pkl for training
                                  and data/multimodaltest-{ontology}.pkl for testing


_____________________________________________________________________________________________________









```

## Contribution

This work currently is under revision in a peer reviewed journal. For use the code or the preprocessed dataset, please open an issue first to discuss what you would like to do. Also you can contact to the corresponding author [Swagarika Jaharlal Giri (swagarika95@gmail.com)](swagarika95@gmail.com)
