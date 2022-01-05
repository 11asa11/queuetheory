from math import log, exp
import random

class Lab_Number_One():
    def __init__(self):
        self.param_lambda = float(input("Введите параметр лямбда: "))
        self.param_n = int(input("Введите параметр n: "))

        self.array_of_x = []
        self.array_of_y = []
        self.array_of_tau = []

        self.__fill_arrays_randomly()
        self.__print_arrays()
        print()

    def __fill_arrays_randomly(self):
        iter = 0
        tau_iter = 0
        for i in range(self.param_n):
            x_iter = random.random()
            y_iter = -(log(1 - x_iter)) / self.param_lambda
            tau_iter += y_iter
            self.array_of_x.insert(iter, x_iter)
            self.array_of_y.insert(iter, y_iter)
            self.array_of_tau.insert(iter, tau_iter)
            iter += 1
        print()

    def __print_arrays(self):
        input_print_arrays = input("Вывести массивы x, y и tau? 1 - Да. 0 - Нет: ")
        if "1" in input_print_arrays:
            count_of_symb = 3
            input_round_num = input("Округлить числа массива до определенного знака после запятой? 3 знака по умолчанию. 1 - Да. 0 - Нет: ")
            if "1" in input_round_num:
                count_of_symb= int(input(
                    "Введите число знаков после запятой: "))

            print("x:  ", end =" ")
            for x in self.array_of_x:
                print(round(x, count_of_symb), end =" ")
            print()

            print("y:  ", end =" ")
            for y in self.array_of_y:
                print(round(y, count_of_symb), end =" ")
            print()

            print("tau:", end =" ")
            for tau in self.array_of_tau:
                print(round(tau, count_of_symb), end =" ")
            print()

    def sort_y(self):
        self.array_of_y.sort()

    def check_through_math_except(self):
        print("Проверка через мат. ожидание: ")
        self.math_except = 1 / self.param_lambda
        self.math_except_new = (1 / self.param_n) * sum(self.array_of_y)
        result_math = abs(self.math_except - self.math_except_new) * 100
        print("Ответ: " + str(round(result_math, 5)) + " %")
        print()

    def check_through_dispersion(self):
        print("Проверка через дисперсию: ")
        disp = 1 / (self.param_lambda * self.param_lambda)
        disp_sum = 0
        for y in self.array_of_y:
            disp_sum += (y - self.math_except_new) * (y - self.math_except_new)
        disp_new = (1 / self.param_n) * disp_sum
        result_disp = abs(disp - disp_new) / disp * 100
        print("Ответ: " + str(round(result_disp, 5)) + " %")
        print()

    def check_through_func_dist(self):
        print("Проверка через функцию распределения: ")
        x_zero = float(input("Введите параметр x0 > 0: "))
        f_x_zero = 1 - exp(-self.param_lambda * x_zero)

        m_x_zero = 0 # Count of intervals which has length smaller than x0
        for length_of_interval in self.array_of_y:
            if (length_of_interval < x_zero):
                m_x_zero += 1

        print("mu(x0) - число интервалов, длина которых меньше x0:  " + str(m_x_zero))
        new_f_x_zero = m_x_zero / self.param_n
        result = abs(f_x_zero - new_f_x_zero)
        print("Ответ: " + str(round(result, 10)))
        print()

def Main():
    print()
    print("Доступные проверки. 1 - через мат. ожидание. 2 - через дисперсию, 3 - через функцию распределения")
    lab = Lab_Number_One()
    exit_flag = True
    while(exit_flag):
        input_checks = input("Введите через запятую необходимые проверки: ")
        if "1" in input_checks:
            lab.check_through_math_except()
        if "2" in input_checks:
            lab.check_through_dispersion()
        if "3" in input_checks:
            lab.check_through_func_dist()

        input_continue = input("Ввести новые параметры? 1 - Да. 0 - Нет: ")
        if "1" in input_continue:
            exit_flag = True
        if "0" in input_continue:
            exit_flag = False

    # SECOND PART
    lab.sort_y()
    m = 1.44 * log(lab.param_n) + 1

    average_length_between_z = int(input("Введите длину разряда: "))
    z_iter = 0
    array_of_z = []
    current_z_value = 0
    for i in range(lab.param_n):
        array_of_z.insert(z_iter, current_z_value)
        z_iter += 1
        current_z_value += average_length_between_z

    array_of_l = []
    for i in range(len(array_of_z) - 1):
        count_of_y_which_length_in_range = 0
        for y_i in lab.array_of_y:
            if y_i >= array_of_z[i] and y_i <= array_of_z[i+1]:
                count_of_y_which_length_in_range += 1
        array_of_l.insert(i, count_of_y_which_length_in_range)

    n = sum(array_of_l)
    array_of_f = []
    for i in range(len(array_of_l)):
        delta_i = array_of_z[i+1] - array_of_z[i]
        array_of_f.insert(i, array_of_l[i] / (n * delta_i))

    print("y: [", end="")
    for y in lab.array_of_y:
        print(round(y, 3), end=", ")
    print()
    print("l: " + str(array_of_l))
    print("z: " + str(array_of_z))
    print("f: " + str(array_of_f))

if __name__ == "__main__":
    Main()