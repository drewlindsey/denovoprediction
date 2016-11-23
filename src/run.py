from pipeline.Pipeline import *

"""
Entry point to the application. Might be temporary?
"""


def main():
    sequence = input("Enter a sequence")
    pipeline = LinearPipeline(sequence)
    pipeline.generate_structure_prediction()


if __name__ == '__main__':
    main()
