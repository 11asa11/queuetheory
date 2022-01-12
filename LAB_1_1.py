from math import log, exp, factorial
import random

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import chi2

def toFixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"

class Lab():
    def __init__(self):
        self.param_lambda = float(input("Введите параметр лямбда: "))
        self.param_n      = int(input("Введите параметр n: "))

        # For 1st part
        self.array_of_x   = []
        self.array_of_y   = []
        self.array_of_tau = []

        # For 2nd part
        self.array_of_z = []
        self.array_of_l = []
        self.array_of_f = []

        # For 4th part
        self.t_zero        = None
        self.m_count_of_intervals = None
        self.k             = None
        self.n_count_table = dict()

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
            count_of_symb = 2
            input_round_num = input("Округлить числа массива до определенного знака после запятой? 2 знака по умолчанию. 1 - Да. 0 - Нет: ")
            if "1" in input_round_num:
                count_of_symb= int(input(
                    "Введите число знаков после запятой: "))

            print("x:  ", end =" ")
            for x in self.array_of_x:
                print(toFixed(x, count_of_symb), end =" ")
            print()

            print("y:  ", end =" ")
            for y in self.array_of_y:
                print(toFixed(y, count_of_symb), end =" ")
            print()

            print("tau:", end =" ")
            for tau in self.array_of_tau:
                print(toFixed(tau, count_of_symb), end =" ")
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

    def clean_arrays_for_second_part(self):
        self.array_of_z.clear()
        self.array_of_l.clear()
        self.array_of_f.clear()

    def first_case_of_bins(self):
        """
        param_m = 1.44 * log(self.param_n) + 1

        average_length_between_z = round(float(input("Введите длину разряда: ")), 2)
        z_iter = 0
        current_z_value = 0
        for i in range(int(param_m)):
            self.array_of_z.insert(z_iter, round(current_z_value, 2))
            z_iter += 1
            current_z_value += round(average_length_between_z, 2)
        """

    def fourth_case_of_bins(self):
        print("Введите диапазоны разрядов, не считая 0")
        self.array_of_z.insert(0, 0.0)
        flag = 1
        z_iter = 1
        while flag == 1:
            flag = int(input("Добавить еще границу разряда? 1 - Да. 0 - Нет: "))
            if flag == 1:
                current_z_value = round(float(input("Введите границу: ")), 2)
                self.array_of_z.insert(z_iter, current_z_value)
                z_iter += 1

    def compute_l_array(self):
        for i in range(len(self.array_of_z) - 1):
            count_of_y_which_length_in_range = 0
            for y_i in self.array_of_y:
                if y_i >= self.array_of_z[i] and y_i < self.array_of_z[i + 1]:
                    count_of_y_which_length_in_range += 1
            self.array_of_l.insert(i, count_of_y_which_length_in_range)

    def compute_f_array(self):
        n = sum(self.array_of_l)
        for i in range(len(self.array_of_l)):
            delta_i = self.array_of_z[i + 1] - self.array_of_z[i]
            self.array_of_f.insert(i, round(self.array_of_l[i] / (n * delta_i), 2))

    def compute_z_average(self):
        z_average = []
        iter = 0
        for i in range(len(self.array_of_z)-1):
            z_average.insert( iter, round((self.array_of_z[i] + self.array_of_z[i+1])/2, 2) )
            iter += 1

        #print("z_average: " + str(z_average))
        return z_average

    def compute_f_ksi_z(self, z_average):
        f_ksi = []
        iter = 0
        for z_i in z_average:
            f_ksi_i = round(self.param_lambda * exp(-self.param_lambda * z_i), 2)
            f_ksi.insert(iter, f_ksi_i)
            iter += 1

        #print("f_ksi: " + str(f_ksi))
        return f_ksi

    def print_arrays_for_second_part(self):
        #print("y: [", end="")

        y_iter = 0
        for y in self.array_of_y:
            if y_iter != (len(self.array_of_y) - 1):
                print(round(y, 2), end=", ")
            else:
                print(round(y, 2),  end="")
            y_iter += 1
        print("]")

        print()
        print("z: " + str(self.array_of_z))
        print("l: " + str(self.array_of_l))
        print("f: " + str(self.array_of_f))

    def compute_mery(self, f_ksi_z):
        mera = float(0)

        for i in range(len(self.array_of_l)):
            mera += (f_ksi_z[i] - self.array_of_f[i])**2

        print("Мера схожести: " + str(round(mera, 5)))

    def create_histogram_first_case(self):
        param_m = int(1.44 * log(self.param_n) + 1)
        #print(param_m)
        bins = plt.hist(self.array_of_y, bins = param_m, density=True, range = (0,max(self.array_of_y)))
        #sns.histplot(data=self.array_of_y, bins=param_m, stat='density')
        iter = 0
        for ranges in bins[1]:
            self.array_of_z.insert(iter, round(ranges, 2))
            iter += 1

        self.compute_l_array()
        self.compute_f_array()
        #self.print_arrays_for_second_part()

        z_average = self.compute_z_average()
        f_ksi_z = self.compute_f_ksi_z(z_average)
        plt.plot(z_average, f_ksi_z, color = "red")

        self.compute_mery(f_ksi_z)

        plt.show()

    def create_histogram_second_case(self):
        param_m = int(1.44 * log(self.param_n) + 1)
        bins_razr = np.interp(np.linspace(0, self.param_n, param_m + 1),
                         np.arange(self.param_n),
                         np.sort(self.array_of_y))
        #print(bins_razr)
        bins = plt.hist(self.array_of_y, bins=bins_razr, density=True)

        iter = 0
        for ranges in bins[1]:
            self.array_of_z.insert(iter, round(ranges, 2))
            iter += 1

        self.compute_l_array()
        self.compute_f_array()

        z_average = self.compute_z_average()
        f_ksi_z = self.compute_f_ksi_z(z_average)
        plt.plot(z_average, f_ksi_z, color = "red")

        self.compute_mery(f_ksi_z)

        plt.show()

    def create_histogram_fourth_case(self):
        plt.hist(self.array_of_y, bins = self.array_of_z, density=True)

        z_average = self.compute_z_average()
        f_ksi_z = self.compute_f_ksi_z(z_average)
        plt.plot(z_average, f_ksi_z, color = "red")

        self.compute_mery(f_ksi_z)

        plt.show()

    def second_part_first_case(self):
        self.clean_arrays_for_second_part()

        self.first_case_of_bins()
        self.create_histogram_first_case()

    def second_part_second_case(self):
        self.clean_arrays_for_second_part()

        self.create_histogram_second_case()

    def second_part_fourth_case(self):
        self.clean_arrays_for_second_part()

        self.fourth_case_of_bins()
        self.compute_l_array()
        self.compute_f_array()
        #self.print_arrays_for_second_part()
        self.create_histogram_fourth_case()

    def third_part(self):
        r = len(self.array_of_l) - 1
        print("Рассматривается гипотеза H0 случайной велечины, распределенной по показательному распределению с")
        print("lamda = " + str(self.param_lambda) + " и r = " + str(r) + " степенями свободы")
        R_0 = 0
        #print(r)
        for i in range(len(self.array_of_l)):
            #print("counter of i: " + str(i) + ". l[i]: " + str(self.array_of_l[i]) + ". z[i]: " + str(self.array_of_z[i]) + ". z[i+1]: " + str(self.array_of_z[i+1]))
            p_i = (1 - exp(-self.param_lambda * self.array_of_z[i + 1])) - (
                        1 - exp(-self.param_lambda * self.array_of_z[i]))
            R_0 += ((self.array_of_l[i] - self.param_n * p_i) ** 2) / (self.param_n * p_i)

        alpha = 1 - round(float(input("Квантиль распределения: ")), 2)
        print("Число степеней свободы = " + str(r))
        critical_value = chi2.ppf(alpha, r)

        if R_0 < critical_value:
            print("R0 = " + str(round(R_0, 2)) + ". Критическое значение = " + str(round(critical_value, 2)) + ". Гипотеза принимается")
        else:
            print("R0 = " + str(round(R_0, 2)) + ". Критическое значение = " + str(round(critical_value, 2)) + ". Гипотеза отвергается")

    def set_t_zero(self):
        upper_limit_t_zero = round(5 / self.param_lambda, 2)
        lower_limit_t_zero = round(3 / self.param_lambda, 2)
        self.t_zero = round(float(input((f'Введите t0 от {lower_limit_t_zero} до {upper_limit_t_zero}: '))),2)

    def compute_intervals(self):
        self.intervals = []
        self.m_count_of_intervals = int(max(self.array_of_tau)/self.t_zero) # is m

        current_value = 0
        for i in range(int(self.m_count_of_intervals + 1)):
            self.intervals.insert(i, round(current_value, 2))
            current_value += self.t_zero

        #print("count_of_intervals: " + str(round(self.m_count_of_intervals, 2)))
        #print("intervals: " + str(intervals))
        return self.intervals

    def count_points(self, intervals):
        array_of_counts = []

        for i in range(len(intervals)-1):
            count = 0
            for tau in self.array_of_tau:
                if tau >= intervals[i] and tau < intervals[i+1]:
                    count += 1
            array_of_counts.insert(i, count)

        #print("array_of_counts: " + str(array_of_counts))

        set_of_counts = set(array_of_counts)
        self.n_count_table = dict.fromkeys(set_of_counts, 0)
        for count in array_of_counts:
            self.n_count_table[count] += 1

        self.k = max(set_of_counts)

        #print("self.k = " + str(self.k))
        #print("set_of_counts: " + str(set_of_counts))
        #print("dict_table: " + str(self.n_count_table))

    def fourth_part_print_table(self):
        print("tau:      ", end=" ")
        for tau in self.array_of_tau:
            print(toFixed(tau, 2), end=" ")
        print()

        print("интервалы:", end=" ")
        for val in self.intervals:
            print(toFixed(val, 2), end=" ")
        print()

        print("[a, b]. a - количество заявок, b - число интервалов")
        for key, value in self.n_count_table.items():
            print("[" + str(key) + "," + str(value) + "]", end=" ")
        print()

    def fourth_part_chi_square(self):
        print("Рассматривается гипотеза H0 случайной велечины, распределенной по закону Пуассона с")
        print("lamda*t0 = " + str(self.param_lambda * self.t_zero) + " и r = " + str(self.k) + " степенями свободы")

        R_0 = 0
        for i in range(self.k + 1):
            #print("counter of i: " + str(i))
            p_i = exp(-self.param_lambda*self.t_zero) * (self.param_lambda * self.t_zero)**i * (1/(factorial(i)))

            n_i = 0
            if i in self.n_count_table:
                n_i = self.n_count_table[i]
            R_0 += ((n_i - self.m_count_of_intervals * p_i) ** 2) / (self.m_count_of_intervals * p_i)

        alpha = 1 - round(float(input("Квантиль распределения: ")), 2)
        print("Число степеней свободы = " + str(self.k) + " = Наибольшее наблюдаемое количество заявок")
        print("Количетсво интервалов (m) = " + str(self.m_count_of_intervals))
        critical_value = chi2.ppf(alpha, self.k)

        if R_0 < critical_value:
            print("R0 = " + str(round(R_0, 2)) + ". Критическое значение = " + str(round(critical_value, 2)) + ". Гипотеза принимается")
        else:
            print("R0 = " + str(round(R_0, 2)) + ". Критическое значение = " + str(round(critical_value, 2)) + ". Гипотеза отвергается")

def snd_part(lab):
    # SECOND PART
    print()
    print("ВТОРАЯ ЧАСТЬ")
    lab.sort_y()
    print("Доступные варианты выбора разрядов:")
    print("1 - Выбор разрядов одинаковой длины")
    print("2 - Выбор разрядов равномерно")
    print("3 - Выбор разрядов в зависимости от точек")
    print("4 - Выбор разрядов вручную")

    exit_flag = True
    while (exit_flag):
        input_choice = input("Введите номер выбора разрядов: ")
        if "1" in input_choice:
            lab.second_part_first_case()
        if "2" in input_choice:
            lab.second_part_second_case()
        if "3" in input_choice:
            lab.check_through_func_dist()
        if "4" in input_choice:
            lab.second_part_fourth_case()

        input_continue = input("Выбрать новый выбор разряда? 1 - Да. 0 - Нет: ")
        if "1" in input_continue:
            exit_flag = True
        if "0" in input_continue:
            exit_flag = False

def fth_part(lab):
    #FOURTH PART
    print()
    print("ЧЕТВЕРТАЯ ЧАСТЬ")
    lab.set_t_zero()
    lab.count_points(lab.compute_intervals())
    print_tab = input("Вывести таблицу количество заявок\число интервалов? 1 - Да. 0 - Нет: ")
    if "1" in print_tab:
        lab.fourth_part_print_table()
    lab.fourth_part_chi_square()

def Main():
    print()
    print("ПЕРВАЯ ЧАСТЬ")
    print("Доступные проверки. 1 - через мат. ожидание. 2 - через дисперсию, 3 - через функцию распределения")
    lab = Lab()
    exit_flag = True

    #FIRST PART
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


    print()
    go_to_part = input("Перейти сразу 4 четвертой части?: 1 - Да. 0 - Нет: ")

    if "1" in go_to_part:
        fth_part(lab)

        snd_part(lab)

        #THIRD PART
        print()
        print("ТРЕТЬЯ ЧАСТЬ")
        lab.third_part()
    else:
        snd_part(lab)

        #THIRD PART
        print()
        print("ТРЕТЬЯ ЧАСТЬ")
        lab.third_part()

        fth_part(lab)

if __name__ == "__main__":
    Main()