from os.path import dirname, abspath
import sys
from multiprocessing import Process, Queue
from gevent.wsgi import WSGIServer

par_path = dirname(dirname(abspath(__file__)))
sys.path.append(par_path)

from flask import Flask, jsonify, render_template, send_from_directory, request, after_this_request
from src.pipeline.Pipeline import *
import threading
import os

project_root = os.path.dirname(__file__)
app = Flask(__name__, template_folder=project_root, static_url_path='/static')
app.debug = True

global pipeline
global thread


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/gen/current')
def get_current_conformation():
    global pipeline
    if pipeline is None:
        return jsonify(result={"status": 400})
    if pipeline.get_current_conformation() is None:
        return jsonify(result={"status": 202})

    pdb = map_conformation_to_pdb(pipeline.get_current_conformation(), app.static_folder, True)

    # pdb = "trythis.pdb"
    print "Sending Conformation Upon Request"
    return send_from_directory(app.static_folder, os.path.basename(pdb))


@app.route('/gen', methods=["POST"])
def gen():
    print "Beginning De Novo Generation"
    data = request.form["sequence"]
    casp_info = get_casp_info(data)
    robetta_dict = casp_info['fragments']
    sequence = casp_info['sequence']
    global pipeline
    global thread
    pipeline = LinearPipeline(data, sequence, robetta_dict)
    # thread = threading.Thread(target=pipeline.generate_structure_prediction(app.static_folder))
    # thread.start()

    run_de_novo()
    return jsonify(result={"status": 200})


@app.route('/gen/complete')
def is_done():
    global pipeline
    if pipeline is None:
        return jsonify(result={"status": 400})

    return jsonify(complete=pipeline.is_complete())


casp_dict = {
    'casp10_t0678': {
        'root': 'CASP10_T0678_UP20725R',
        'sequence': '',
        'fragments': ''
    },
    'casp10_t0651': {
        'root': 'CASP10_T0651_LGR82',
        'sequence': '',
        'fragments': ''
    },
    'casp10_t0694': {
        'root': 'CASP10_T0694_APC100075',
        'sequence': '',
        'fragments': ''
    },
    'casp10_t0757': {
        'root': 'CASP10_T0757_APC103790',
        'sequence': '',
        'fragments': ''
    },
    'casp10_t0666': {
        'root': 'CASP10_T0666_UCI_BBCS',
        'sequence': '',
        'fragments': ''
    },
    'casp11_t0837': {
        'root': 'CASP11_T0837_YPO2654',
        'sequence': '',
        'fragments': ''
    },
    'casp11_t0792': {
        'root': 'CASP11_T0792_Oskar-N',
        'sequence': '',
        'fragments': ''
    },
    'casp11_t0806': {
        'root': 'CASP11_T0806_YaaA',
        'sequence': '',
        'fragments': ''
    },
    'casp11_t0843': {
        'root': 'CASP11_T0843_Ats13',
        'sequence': '',
        'fragments': ''
    },
    'casp11_t0856': {
        'root': 'CASP11_T0856_HERC1',
        'sequence': '',
        'fragments': ''
    }
}


# def generate_pdb():
#    global pipeline
#    global curr_pdb
#    if


def run_de_novo():
    global pipeline
    print "%%%%% STARTING %%%%%%"
    gen_process = Process(target=pipeline.generate_structure_prediction(app.static_folder))
    gen_process.start()
    print "%%%%% HOPEFULLY THIS SHOWS UP NEXT %%%%%%"
    # pdb_process = Process(target=generate_pdb(queue))
    # pdb_process.start()


def get_casp_info(casp_name):
    seq = ''
    path3mer = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 'aat000_03_05.200_v1_3.txt')
    path9mer = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 'aat000_09_05.200_v1_3.txt')
    path_fasta = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 't000_.fasta')
    with open(path_fasta) as fasta:
        next(fasta)
        for line in fasta:
            seq += line[:-1]  # remove newline char
    return {
        'sequence': seq,
        'fragments': {
            3: path3mer,
            9: path9mer
        }
    }


if __name__ == "__main__":
    http_server = WSGIServer(('', 5000), app)
    http_server.serve_forever()
    # app.run(debug=True, threaded=True)
