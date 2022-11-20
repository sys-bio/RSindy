
import reactive_sindy
import sys

if __name__ == "__main__":

    # todo: separate data generation

    # sys.argv[1]: model name
    # sys.argv[2]: end time (start time is always 0)
    # sys.argv[3]: time points

    output = reactive_sindy.rsindy(sys.argv[1], sys.argv[2], sys.argv[3])

    print()
    print('number of iterations:', output[0])
    print()
    print('alpha parameter:', output[1])
    print()
    print('ansatz reactions')
    print(output[2])
    print()
    print('inferred parameter values')
    for each in output[3]:
        print(each)
    print()
