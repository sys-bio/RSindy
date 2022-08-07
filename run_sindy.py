
import reactive_sindy

if __name__ == "__main__":

    # import test model and run the algorithm
    # todo: separate data generation
    # todo: allow constraints on the ansatz reactions
    output = reactive_sindy.rsindy('lv.txt')

    print()
    print('number of iterations:', output[0])
    print()
    print('alpha parameter:', output[1])
    print()
    print('ansatz reactions')
    print(output[2])
    # print()
    print('inferred parameter values')
    for each in output[3]:
        print(each)
    print()
