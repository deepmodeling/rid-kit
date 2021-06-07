from google.protobuf import text_format
from tensorflow.python.platform import gfile
import tensorflow as tf


def convert_pb_to_pbtxt(pbfile, pbtxtfile):
    with gfile.FastGFile(pbfile, 'rb') as f:
        graph_def = tf.compat.v1.GraphDef()
        graph_def.ParseFromString(f.read())
        tf.import_graph_def(graph_def, name='')
        tf.io.write_graph(graph_def, './', pbtxtfile, as_text=True)


def convert_pbtxt_to_pb(pbtxtfile, pbfile):
    with tf.gfile.FastGFile(pbtxtfile, 'r') as f:
        graph_def = tf.compat.v1.GraphDef()
        file_content = f.read()
        # Merges the human-readable string in `file_content` into `graph_def`.
        text_format.Merge(file_content, graph_def)
        tf.io.write_graph(graph_def, './', pbfile, as_text=False)


if __name__ == "__main__":
    convert_pb_to_pbtxt("./benchmark_case/000/graph.001.pb",
                        "./benchmark_case/000/graph.001.pbtxt")
