from flask import Flask
from flask import jsonify,request
from flask_cors import CORS, cross_origin
app = Flask(__name__)
CORS(app)
import sys,os
import DataGenerationAPI 



@app.route('/accession')
@cross_origin()
def get_accession_number_status():
    if 'accession_no' in request.args and 'ontology' in request.args:
        try:
            data=DataGenerationAPI.analyze_accession_status(request.args['accession_no'],request.args['ontology'])
            message = {
                'status': 200,
                    'message': 'OK',
                    'data': data
                    }
            resp = jsonify(message)
            resp.status_code = 200
            print(resp)
            return resp
        except:
            message = {
                    'status': 500,
                        'message': 'Error has occured',
                        'data': {}
                        }
            resp = jsonify(message)
            print(resp)
            return resp

    else:
        return 'Sorry Accession Number is not received'







if __name__ == "__main__":
    app.run(debug = True,host = "172.16.26.35", port = 5003)