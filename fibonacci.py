'''
Fibonacci
Write a program that prints the Fibonacci numbers from n to m (inclusive).

Requirements

for each value print corresponding fibonacci number
1 <= n < m <= 250
optimal time complexity

Input
Two numbers in two lines (n, m)

Output
One result per line according to requirements

Example
Sample input:
20
25
Sample output:
6765
10946
17711
28657
46368
75025
'''

def return_fib_num(v):
    if v == 0:
        return 0
    elif v == 1:
        return 1
    else:
        return return_fib_num(v-1) + return_fib_num(v-2)

def main():
    ### obtain inputs
    try:
        num1 = int(input())
        num2 = int(input())
    except ValueError:
        print('Error1: Please provide integers between 1 and 250 as inputs. The second integer must be larger than the first integer.')
        raise

    ### check input requirements
    try:
        flag = False
        if num1 >= 1 and num2 <= 250 and num1 < num2:
            flag = True
        assert flag
    except AssertionError:
        print('Error2: Please provide integers between 1 and 250 as inputs. The second integer must be larger than the first integer.')
        raise

    ### print the fibonacci numbers
    for num in range(num1, num2+1):
        print(return_fib_num(num))

    return

if __name__ == '__main__':
    main()

