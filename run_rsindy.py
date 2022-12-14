
import reactive_sindy
import sys
import time

if __name__ == "__main__":

    # todo: separate data generation

    # Note: the model file, ansatz antimony file, and ansatz python module must all be in the same directory
    # as the run_sindy and reactive_sindy scripts

    # sys.argv[1]: model name
    # sys.argv[2]: simulation end time (start time is always 0)
    # sys.argv[3]: number of simulation points

    initial = time.time()
    output = reactive_sindy.rsindy(sys.argv[1], sys.argv[2], sys.argv[3])
    final = time.time()

    print()
    print('time: ', final - initial)
    print()
    print('number of iterations:', output[0])
    print()
    print('alpha parameter:', output[1])
    print()
    print('inferred parameter values')
    for each in output[3]:
        print(each)
    print()
