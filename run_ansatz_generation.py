
import ansatz_generation
import sys

if __name__ == "__main__":

    # sys.argv[1]: model name
    # sys.argv[2]: number of ansatz reactions
    # sys.argv[3]: number of ansatz sets

    # This can subset the full set of ansatz reactions

    ansatz_generation.generate_ansatz(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
