from os.path import dirname, abspath
import sys
from multiprocessing import Process, Queue

import time
from gevent.wsgi import WSGIServer
from celery import Celery

par_path = dirname(dirname(abspath(__file__)))
sys.path.append(par_path)

from flask import Flask, jsonify, render_template, send_from_directory, request, url_for
from src.pipeline.Pipeline import *
import threading
import os

project_root = os.path.dirname(__file__)
app = Flask(__name__, template_folder=project_root, static_url_path='/static')
app.debug = True

# Celery configuration
app.config["CELERY_BROKER_URL"] = "redis://localhost:6379/0"
app.config["CELERY_RESULT_BACKEND"] = "redis://localhost:6379/0"

# Initialize Celery
celery = Celery(app.name, broker=app.config["CELERY_BROKER_URL"])
celery.conf.update(app.config)

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
    # print "Sending Conformation Upon Request"
    return send_from_directory(app.static_folder, os.path.basename(pdb))


@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = generate_conformation.AsyncResult(task_id)
    # print task.state
    #if task.state == "NEXT":
     #   return jsonify({"current": task.info["current"],
     #                   "total": task.info["total"]})
    if task.state == "PDB_CHANGE":
        return send_from_directory(app.static_folder, os.path.basename(task.info["minimum"]))

    return jsonify({"status": task.state})


@app.route('/generate', methods=["POST"])
def generate():
    name = request.form["sequence"]
    casp_info = get_casp_info(name)
    robetta_dict = casp_info["fragments"]
    sequence = casp_info["sequence"]
    experimental = casp_info["experimental_pdb"]
    experimental2 = casp_info["experimental_pdb2"]

    task = generate_conformation.apply_async((name, robetta_dict, sequence, experimental2))

    return jsonify({}), 202, {'Location': url_for('taskstatus', task_id=task.id)}


@celery.task(bind=True)
def generate_conformation(self, name, robetta_dict, sequence, experimental):
    """Background task that runs to generate the Conformation with frequent
    updates in the form of PDB files and other data"""
    global pipeline
    # pipeline = LinearPipeline(name, sequence, robetta_dict)
    new_dict = {}
    for key in robetta_dict:
        new_dict[int(key)] = robetta_dict[key]
    frag_lib = RobettaFragmentLibrary(sequence)
    frag_lib.generate(new_dict)
    seef = TMScore()#DFirePotential()
    score = {"rmsd": RMSDScore(), "tm-score": TMScore()}
    self.conformation = LinearBackboneConformation(name, sequence, experimental)
    self.conformation.initialize()
    sampler = ConformationSampler(self.conformation, experimental, seef, frag_lib, app.static_folder, score)

    count = 0

    while sampler.has_next():
        old_pdb = self.conformation.get_pdb_file()
        self.conformation = sampler.next_conformation()
        new_pdb = self.conformation.get_pdb_file()

        self.update_state(state="NEXT",
                          meta={"current": count,
                                "total": sampler.get_k_max()})
        count += 1
        if old_pdb != new_pdb:
            #print "New PDB"
            curr_en = sampler.get_current_energy()
            print "CURRENT TMSCORE: " + str(curr_en)
            curr_score = sampler.get_current_score()
            print "CURRENT RMSD: " + str(curr_score)
            best_en = sampler.get_best_energy()
            print "BEST TMSCORE: " + str(best_en)
            score_for_best_tm = sampler.get_best_score_for_tm()
            print "RMSD FOR TMSCORE: " + str(score_for_best_tm)
            best_rmsd = sampler.get_best_score()
            print "BEST RMSD: " + str(best_rmsd)
            tm_for_best_score = sampler.get_best_tm_for_score()
            print "TMSCORE FOR RMSD: " + str(tm_for_best_score)
            
            self.update_state(state="PDB_CHANGE",
                              meta={"minimum": sampler.minimum().get_pdb_file()})

            # TODO submit to 3dmol.js

    best_dir = os.path.join(app.static_folder, "best")
    if not os.path.exists(best_dir):
        os.makedirs(best_dir)

    timestr = time.strftime("%Y_%m_%d-%H_%M_%S")
    file_name = os.path.join(best_dir, name + "_" + timestr + ".pdb")
    with open(file_name, 'w+') as min_file:
        with open(sampler.minimum().get_pdb_file()) as pdb_file:
            for line in pdb_file:
                min_file.write(line)

    self.result = sampler.score_conformation()
    self.conformation = sampler.minimum()
    self.update_state(state="PDB_CHANGE",
                      meta={"minimum": sampler.minimum().get_pdb_file()})

    print "%%%%%%%%%%%%%%%%%%%%% FINAL %%%%%%%%%%%%%%%%%%%%%"
    print self.conformation.get_pdb_file()
    for key, value in self.result:
        print "%%%%%%%%%%%%%%% [" + str(key) + "] " +  str(self.result) + " %%%%%%%%%%%%%%%"
    
    # print "Beginning while loop"
    # while pipeline.is_complete():
    #     if random.random() < 0.10:
    #         self.update_state(state="PROGRESS",
    #                           meta={'complete': pipeline.is_complete(),
    #                                 'pdb':  pipeline.get_current_conformation().get_pdb_file()})
    #          time.sleep(1)
    return {"status": 200, "result": self.result}


"""@app.route('/gen', methods=["POST"])
def gen():
    print "Beginning De Novo Generation"
    data = request.form["sequence"]
    casp_info = get_casp_info(data)
    robetta_dict = casp_info['fragments']
    sequence = casp_info['sequence']
    experimental_pdb = casp_info['experimental_pdb']
    global pipeline
    global thread
    pipeline = LinearPipeline(data, sequence, robetta_dict, experimental_pdb)
    # thread = threading.Thread(target=pipeline.generate_structure_prediction(app.static_folder))
    # thread.start()

    # run_de_novo()
    return jsonify(result={"status": 200})


@app.route('/gen/complete')
def is_done():
    global pipeline
    if pipeline is None:
        return jsonify(result={"status": 400})

    return jsonify(complete=pipeline.is_complete())"""

casp_dict = {
    'casp10_t0678': {
        'root': 'CASP10_T0678_UP20725R',
        'sequence': '',
        'fragments': '',
        'expPDB': '4epz.pdb', 
        'experimental_pdb': '',
        'casp': 'CASP10',
        'expPDB2': 'T0678.pdb'
    },
    'casp10_t0651': {
        'root': 'CASP10_T0651_LgR82',
        'sequence': '',
        'fragments': '',
        'expPDB': '4f67.pdb',
        'experimental_pdb': '',
        'casp': 'CASP10',
        'expPDB2': 'T0651.pdb'
    },
    'casp10_t0694': {
        'root': 'CASP10_T0694_APC100075',
        'sequence': '',
        'fragments': '',
        'expPDB': '5jh8.pdb',
        'experimental_pdb': '',
        'casp': 'CASP10',
        'expPDB2': 'T0694.pdb'
    },
    'casp10_t0757': {
        'root': 'CASP10_T0757_APC103790',
        'sequence': '',
        'fragments': '',
        'expPDB': '4gak.pdb',
        'experimental_pdb': '',
        'casp': 'CASP10',
        'expPDB2': 'T0757.pdb'
    },
    'casp10_t0666': {
        'root': 'CASP10_T0666_UCI_BBCS',
        'sequence': '',
        'fragments': '',
        'expPDB': '3ux4.pdb',
        'experimental_pdb': '',
        'casp': 'CASP10',
        'expPDB2': 'T0666.pdb'
    },
    'casp11_t0837': {
        'root': 'CASP11_T0837_YPO2654',
        'sequence': '',
        'fragments': '',
        'expPDB': '5tf3.pdb',
        'experimental_pdb': '',
        'casp': 'CASP11',
        'expPDB2': 'T0837.pdb'
    },
    'casp11_t0792': {
        'root': 'CASP11_T0792_Oskar-N',
        'sequence': '',
        'fragments': '',
        'expPDB': '5a49.pdb',
        'experimental_pdb': '',
        'casp': 'CASP11',
        'expPDB2': 'T0792.pdb'
    },
    'casp11_t0806': {
        'root': 'CASP11_T0806_YaaA',
        'sequence': '',
        'fragments': '',
        'expPDB': '5caj.pdb',
        'experimental_pdb': '',
        'casp': 'CASP11',
        'expPDB2': 'T0806.pdb'
    },
    'casp11_t0843': {
        'root': 'CASP11_T0843_Ats13',
        'sequence': '',
        'fragments': '',
        'expPDB': '4xau.pdb',
        'experimental_pdb': '',
        'casp': 'CASP11',
        'expPDB2': 'T0843.pdb'
    },
    'casp11_t0856': {
        'root': 'CASP11_T0856_HERC1',
        'sequence': '',
        'fragments': '',
        'expPDB': '4qt6.pdb',
        'experimental_pdb': '',
        'casp': 'CASP11',
        'expPDB2': 'T0856.pdb'
    }
}


def get_casp_info(casp_name):
    seq = ''
    path3mer = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 'aat000_03_05.200_v1_3.txt')
    path9mer = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 'aat000_09_05.200_v1_3.txt')
    path_fasta = os.path.join(app.static_folder, casp_dict[casp_name]['root'], 't000_.fasta')
    path_exp = os.path.join(app.static_folder, 'ExperimentalNativeStructures', casp_dict[casp_name]['root'], casp_dict[casp_name]['expPDB'])
    path_exp2 = os.path.join(app.static_folder, 'NativeStructures', casp_dict[casp_name]['casp'], casp_dict[casp_name]['expPDB2'])
    with open(path_fasta) as fasta:
        next(fasta)
        for line in fasta:
            seq += line[:-1]  # remove newline char
    return {
        'sequence': seq,
        'fragments': {
            3: path3mer,
            9: path9mer
        },
        'experimental_pdb': path_exp,
        'experimental_pdb2': path_exp2
    }
    
if __name__ == "__main__":
    http_server = WSGIServer(('', 5000), app)
    http_server.serve_forever()
    # app.run(debug=True, threaded=True)
