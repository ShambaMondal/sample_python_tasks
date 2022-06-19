'''
FizzBuzz
Write a program that prints the integers from n to m (inclusive).
Requirements:

- for multiples of three, print Fizz (instead of the number)
- for multiples of five, print Buzz (instead of the number)
- for multiples of both three and five, print FizzBuzz (instead of the number)
and 1 <= n < m <= 10000


Input
Two numbers in two lines (n, m)

Output
One result per line including requirements

Example
Sample input:
7
16
Sample output:
7
8
Fizz
Buzz
11
Fizz
13
14
FizzBuzz
16
'''
### obtain inputs
def obtain_inputs():
    try:
        num1 = int(input())
        num2 = int(input())
    except ValueError:
        print('Error1: Please provide integers between 1 and 10000 as inputs. The second integer must be larger than the first integer.')
        raise

    ### check input requirements
    try:
        flag = False
        if num1 >= 1 and num2 <= 10000 and num1 < num2:
            flag = True
        assert flag
    except AssertionError:
        print('Error2: Please provide integers between 1 and 10000 as inputs. The second integer must be larger than the first integer.')
        raise
    return num1, num2

### print
def print_output(num1, num2):
    for num in range(num1, num2+1):
        mult3 = (num%3 == 0)
        mult5 = (num%5 == 0)

        value_dict = {(True, True): 'FizzBuzz', (True, False): 'Fizz', (False, True): 'Buzz', (False, False): num}
        print(value_dict[(mult3, mult5)])
    return

def main():
    n1, n2 = obtain_inputs()
    print_output(n1, n2)

if __name__ == '__main__':
    main()
