
import reactive_sindy
import sys
import time

if __name__ == "__main__":

    # todo: separate data generation

    # sys.argv[1]: model name
    # sys.argv[2]: simulation end time (start time is always 0)
    # sys.argv[3]: number of simulation points

    initial = time.time()
    output = reactive_sindy.rsindy('lv.txt', 15, 200)
    final = time.time()

    print()
    print('time:', final - initial)
    print()
    print('number of iterations:', output[0])
    print()
    print('alpha parameter:', output[1])
    print()
    print('inferred parameter values')
    for each in output[3]:
        print(each)
    print()
