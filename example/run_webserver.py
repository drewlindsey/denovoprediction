from flask import Flask, jsonify, render_template, send_from_directory
from src.pipeline.Pipeline import *
import threading
import os

project_root = os.path.dirname(__file__)
app = Flask(__name__, template_folder=project_root, static_url_path='/')


global pipeline
global thread


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/gen/current')
def get_current_conformation():
    print("test3")
    global pipeline
    # pdb = pipeline.get_pdb_file()
    pdb = "trythis.pdb"
    return send_from_directory(app.root_path, pdb)


@app.route('/gen')
def gen(sequence):
    print("test2")
    global pipeline
    global thread
    pipeline = LinearPipeline(sequence)
    thread = threading.Thread(target=pipeline.generate_structure_prediction())


@app.route('/gen/continue')
def is_done():
    global pipeline
    print("test")
    # return jsonify(complete=pipeline.is_complete());
    return jsonify(complete=False)


if __name__ == "__main__":
    app.run()
