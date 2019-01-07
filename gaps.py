#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 4 17:48:32 2018

@author: saadiq 
"""
from fractions import Fraction 
import math
import sys
import os
from pulp import *

class Fr:
    def __init__(self, n = 1, d = 1):
        self.numerator = n
        self.denominator = d
        
    def get_fraction(self):
        return Fraction(self.numerator, self.denominator)
    
    def __repr__(self):
        return '{}/{}'.format(self.numerator, self.denominator)

class Term:
    def __init__(self, l, new_var):
        self.coeff = l
        self.variable = new_var

    def __repr__(self):
        return ""+"("+str(self.coeff)+")"+str(self.variable)

    def equals(self, other):
        return self.coeff == other.coeff and self.variable == other.variable

class Equation:
    def __init__(self, new_arg1, new_arg2):
        self.arg1 = new_arg1
        self.arg2 = new_arg2

    def __repr__(self):
        string = ""
        if len(self.arg1) != 0:
            for i in range(len(self.arg1) - 1):
                string += str(self.arg1[i])
                string += " + "
            string += str(self.arg1[len(self.arg1) - 1])
        string += " = "
        if len(self.arg2) != 0:
            for i in range(len(self.arg2) - 1):
                string += str(self.arg2[i])
                string += " + "
            string += str(self.arg2[len(self.arg2) - 1])

        return string

class ThreeShareException(Exception):
    pass

class BillCaseException(Exception):
    pass


class Interval: 
    type_list = [0,0,0]
      
    def __init__(self, new_num1=0, new_num2=0, new_type=0, new_num_shares=0):
        self.low = min(new_num1, new_num2)
        self.high = max(new_num1, new_num2)
        self.type = new_type
        self.num_shares = new_num_shares
        
    def contains_bound(self, bound):
        return bound >= self.low and bound <= self.high
    
    def has_bound(self, bound):
        return bound == self.low or bound == self.high
    
    def equals(self, other):
        return other.low == self.low and other.high == self.high
    
    def __repr__(self):
        if self.type == 1 or self.type == 2:
            return "(%s, %s)" % (self.low, self.high)
        else:
            return "(%s [no shares] %s)" % (self.low, self.high)

                
class Student:
    def __init__(self, nums, share_type):
        self.numbers = nums
        self.share_type = share_type
        
    def equals(self, s):
        return self.numbers == s.numbers 
    
    def contains(self, int_type):
        inside = False
        for elem in self.numbers:
            if elem == int_type:
                inside = True
                
        return inside
    
    def num_occurences(self, int_type):
        count = 0
        for elem in self.numbers:
            if elem == int_type:
                count += 1
        
        return count
    
    def __repr__(self):
        new_str = "("
        nums = list(map(lambda x: str(x), self.numbers))
        new_str += " ".join(nums)
        new_str += ")"
        
        return new_str

class IntervalManager:
    def __init__(self):
        self.interval = []
        
    def find_interval(self, itvl):
        index = -1
        for i in range(len(self.interval)):
            if self.interval[i].equals(itvl):
                index = i
        
        return index
    
    def find_bound(self, bound):
        index = -1
        for i in range(len(self.interval)):
            if self.interval[i].contains_bound(bound):
                index = i
                
        return index
    
    def add_ignore_everything(self, index, itvl):
        self.interval.insert(index, itvl)
        
    def check_same_bounds(self):
        counter = 0
        while counter < len(self.interval):
            if self.interval[counter].low == self.interval[counter].high:
                del self.interval[counter]
            else:
                counter += 1
                
    def get_highest_bound(self):
        return self.interval[len(self.interval) - 1].high
    
    def get_lowest_bound(self):
        return self.interval[0].low
    
    def length(self):
        return len(self.interval)
    
    def get(self, index):
        return self.interval[index]
    
    def int_print(self):
        for elem in self.interval:
            print(elem, " ", end='')
        print()
        
    def get_all(self): 
        return self.interval
    
    def get_all_not(self, bad):
        new_list = []
        for elem in self.interval:
            if elem.type != bad:
                new_list.append(elem)
                
        return new_list
    
    def find_buddy(self, itvl, denom):
        equal = Interval(denom - itvl.high, denom - itvl.low, itvl.type)
        index = -1
        for i in range(len(self.interval)):
            if self.interval[i].equals(equal): 
                index = i
        
        return index
    
    def get_all_type(self, index_type): #getAll(double type) in Java code
        lst = []
        
        for elem in self.interval:
            if elem.type == index_type:
                lst.append(elem)
                
        return lst
        
    def remove(self, index):
        del self.interval[index]
                
    def add_ignore_share(self, itvl):
        low_interval_index = self.find_bound(itvl.low)
        high_interval_index = self.find_bound(itvl.high)
        
        if low_interval_index == high_interval_index == -1:
            self.interval.append(itvl)
        elif len(self.interval) == 1 and itvl.low < self.interval[0].high and itvl.low > self.interval[0].low and self.interval[0].type != 0: 
            self.interval.append(itvl)
        elif len(self.interval) == 1 and itvl.high > self.interval[0].low and itvl.high < self.interval[0].high and self.interval[0].type != 0:
            self.interval.insert(0, itvl)
        elif high_interval_index == low_interval_index and len(self.interval) == 1 and (itvl.high > self.interval[0].high or itvl.low < self.interval[0].low or (high_interval_index == -1 or low_interval_index == -1) and len(self.interval) == 1 and (itvl.high > self.interval[0].high or itvl.low < self.interval[0].low)):
            if itvl.high > self.interval[0].high:
                if itvl.low < self.interval[0].low:
                    del self.interval[0]
                    self.interval.append(Interval(itvl.low, itvl.high, itvl.type))
                else:
                    self.interval[0].high = itvl.low
                    self.interval.append(itvl)
            else:
                if itvl.high > self.interval[0].high:
                    del self.interval[0]
                    self.interval.append(Interval(itvl.low, itvl.high, itvl.type))
                else:
                    self.interval[0].low = itvl.high
                    self.interval.insert(0, itvl)
            self.check_same_bounds()
        elif high_interval_index == low_interval_index:   
            if self.interval[low_interval_index].type != 0 and itvl.type != 0:
                itvl.type = self.interval[low_interval_index].type
            lowItvl = Interval(self.interval[low_interval_index].low, itvl.low, self.interval[low_interval_index].type)
            highItvl = Interval(itvl.high, self.interval[high_interval_index].high, self.interval[high_interval_index].type) 
            
            del self.interval[high_interval_index]
            self.interval.insert(high_interval_index, lowItvl)

            self.interval.insert(high_interval_index + 1, itvl)
            
            self.interval.insert(high_interval_index + 2, highItvl)
            self.check_same_bounds()
        else:
            if high_interval_index - low_interval_index > 1:
                num_iterations = high_interval_index - low_interval_index - 1
                counter = 1
                start = low_interval_index + 1
                while counter <= num_iterations:
                    del self.interval[start]
                    counter += 1
                
                self.interval.insert(low_interval_index + 1, itvl)
                self.interval[low_interval_index].high = itvl.low
                self.interval[low_interval_index + 2].low = itvl.high
                self.check_same_bounds()
            else:
                self.interval[low_interval_index].high = itvl.low 
                self.interval[high_interval_index].low = itvl.high
                self.interval.insert(low_interval_index + 1, itvl)
                self.check_same_bounds()
                   
    
    def add_interval(self, itvl):
        low_interval_index = self.find_bound(itvl.low)
        high_interval_index = self.find_bound(itvl.high)
        
        if (low_interval_index == -1 or high_interval_index == -1):
            self.interval.append(itvl)
        elif (low_interval_index == high_interval_index): 
            if (self.interval[low_interval_index].type != 0 and itvl.type != 0):
                itvl.type = self.interval[low_interval_index].type
            low_itvl = Interval(self.interval[low_interval_index].low, itvl.low, self.interval[low_interval_index].type)
            high_itvl = Interval(itvl.high, self.interval[high_interval_index].high, self.interval[high_interval_index].type)
            
            del self.interval[high_interval_index]
            self.interval.insert(high_interval_index, low_itvl)
            
            self.interval.insert(high_interval_index + 1, itvl)
            
            self.interval.insert(high_interval_index + 2, high_itvl)
            self.check_same_bounds()
        else:
            if (high_interval_index - low_interval_index > 1):
                the_type = 0
                for i in range(low_interval_index, high_interval_index + 1):
                    if (self.interval[i].type != 0):
                        the_type = self.interval[i].type
                num_iterations = high_interval_index - low_interval_index - 1
                counter = 1
                start_index = low_interval_index + 1
                while (counter <= num_iterations):
                    del self.interval[start_index]
                    counter += 1
                if itvl.type != 0:
                    itvl.type = the_type
                self.interval.insert(low_interval_index + 1, itvl)
                self.interval[low_interval_index].high = itvl.low
                self.interval[low_interval_index + 2].low = itvl.high;
                self.check_same_bounds()
            else:
                the_type = 0
                for i in range(low_interval_index, high_interval_index + 1):
                    if self.interval[i].type != 0:
                        the_type = self.interval[i].type
                if itvl.type != 0:
                    itvl.type = the_type
                    
                self.interval[low_interval_index].high = itvl.low
                self.interval[high_interval_index].low = itvl.high
                self.interval.insert(low_interval_index + 1, itvl)
                self.check_same_bounds()
            
class MuffinsAuto:
    def __init__(self, m, s, alpha_num, alpha_den):
        self.m = m
        self.s = s
        
        self.is_more_gaps = True
        self.share = Fr(m, s)
        self.alpha = Fr(alpha_num, alpha_den)
        self.manager = IntervalManager()
        
        self.low_share = 1
        self.high_share = 1
        self.num_low_share = 1
        self.num_high_share = 1
        self.turn = 1
        self.low = 0
        self.high = 0

        self.easy_erik = 0
        
        if alpha_den % self.share.denominator == 0:
            self.share = self.convert(self.share, alpha_den) 
        else:
            self.lcm = self.lcm(alpha_den, s) 
            self.share = self.convert(self.share, self.lcm)
            self.alpha = self.convert(self.alpha, self.lcm)
            
        self.calcS()
        
    def convert(self, a, common):
        factor = common//a.denominator
        res = Fr(a.numerator * factor, common)

        return res

    def lcm(self, a, b):
        return (a*b)//math.gcd(a, b)   

    
    def combination_util(self, arr, data, start, end, index, r, students):
        if index == r:
            for j in range(r):
                students.append(data[j])
            return
        
        i = start
        while (i <= end and end - i + 1 >= r - index):
            data[index] = arr[i]
            self.combination_util(arr, data, i+1, end, index+1, r, students)
            i += 1
            
    
    def print_combination(self, arr, n, r, students):
        data = [0]*r
        
        self.combination_util(arr, data, 0, n-1, 0, r, students)
        
    
    def combination(self, num_intervals, share_type):
        students = []
        
        arr = [0]*(num_intervals*share_type)
        
        interval_count = 1
        for i in range(len(arr)):
            arr[i] = interval_count
            if (i+1) % share_type == 0:
                interval_count += 1
                
        self.print_combination(arr, len(arr), share_type, students)
        return students
        
    
    def three_share_test(self):
        super_high = self.high_share + 1
        super_low = self.low_share - 1
        super_high_fraction = Fr(self.m, self.s * super_high)
        super_low_fraction = Fr(self.m, self.s * super_low)
        
        is_less_one = (super_high_fraction.get_fraction() - self.alpha.get_fraction()).numerator > 0
        is_less_two = (Fraction(1) - super_low_fraction.get_fraction() - self.alpha.get_fraction()).numerator > 0               
        
        if is_less_one or is_less_two:
            raise ThreeShareException("DARN, THIS IS THE CASE BILL WARNED ME ABOUT WHERE THERE ARE THREE NUMBERS-OF-SHARES.")
       
         
    def calcS(self):
        pieces = 2 * self.m
        self.low_share = pieces//self.s
        self.high_share = self.low_share + 1
        
        num1 = pieces - self.high_share*self.s
        num2 = self.low_share*self.s - pieces
        denom = self.low_share - self.high_share
        
        x = num1//denom
        y = num2//denom
        
        if (denom == 0):
            print("\nThe equation has no solution.")
        else:
            self.num_low_share = x
            self.num_high_share = y
            
        alpha_decimal = self.alpha.numerator/self.alpha.denominator
        
        print("Try to prove f(%s, %s) <= %s (%s)" % (self.m, self.s, self.alpha, alpha_decimal))
        
        if (self.m != self.share.numerator):
            print("Note that %s/%s = %s" % (self.m, self.s, self.share))    
        print("s%s = %s" %(self.low_share, self.num_low_share))
        print("s%s = %s" %(self.high_share, self.num_high_share))
        
        print("All numbers are assumed to have a denominator of %s" % (self.alpha.denominator))
        
        Interval.type_list = [0, self.low_share, self.high_share]
        
    def subtract_from(self, a, b):
        if a.denominator == 0 or b.denominator == 0:
            raise ValueError("Invalid denominator")
        
        common = a.denominator
        
        commonA = self.convert(a, common)
        commonB = self.convert(b, common)
        
        diff = Fr(commonB.numerator - commonA.numerator, common)
        
        return diff
        
    def find_upper_range(self):
        new_bound = (-1 *((self.low_share - 1)*self.manager.get_highest_bound() - self.share.numerator))
        
        if new_bound < self.low:
            new_bound = self.low
            
        res = Interval(new_bound, self.subtract_from(self.alpha, Fraction(1)).numerator, 1, self.num_low_share * self.low_share)
        return res

    
    def find_lower_range(self):
        new_bound = (-1 *((self.high_share - 1)*self.manager.get_lowest_bound() - self.share.numerator))
        
        if new_bound > self.high:
            new_bound = self.high
            
        self.manager.add_ignore_share(Interval(self.alpha.numerator, new_bound, 2, self.num_high_share * self.high_share))
    
        return new_bound
    
    
    def set_interval(self):
        self.low = self.alpha.numerator
        res = self.subtract_from(self.alpha, Fraction(1))
        self.high = res.numerator

        self.manager.add_interval(Interval(self.low, self.high, 0))
        
    
    def find_blocks(self):
        temp_list = []
        
        for i in range(self.manager.length()):
            if (self.manager.get(i).type == 0):
                temp_list.append(self.manager.get(i))
                
        for elem in temp_list:
            self.manager.add_interval(Interval(self.alpha.denominator - elem.high, self.alpha.denominator - elem.low, 0))
    
    
    def myint(self, x):
        if x != int(x):
            return x
        else:
            return int(x) 

    
    def mid_point(self):
        mid = self.myint((self.manager.get_lowest_bound() + self.manager.get_highest_bound())/2) 
        index = self.manager.find_bound(mid)
        elem = self.manager.get(index)    
        itvl = Interval(elem.low, mid, elem.type, elem.num_shares)    
        self.manager.add_interval(itvl)
    
        
    def test(self, possible, share_type): 
        lst = []
        temp_list = []
        interval_list = []
        
        if share_type == self.high_share:
            interval_list = self.manager.get_all_type(2)
        else:
            interval_list = self.manager.get_all_type(1)
            
        i = 0
        while (i < len(possible) - share_type + 1):
            total_max = 0
            total_min = 0
            
            for j in range(i, i + share_type):
                temp_list.append(possible[j])
                total_max += interval_list[possible[j] - 1].high
                total_min += interval_list[possible[j] - 1].low
                
            if total_max > self.share.numerator and total_min < self.share.numerator:
                for elem in temp_list:
                    lst.append(elem)
                    
            temp_list = []
            i += share_type
                    
        return lst
        
    
    def find_duplicates(self, list1):
        final_list = []
        
        for i in range(len(list1)):
            has_equals = False
            
            for elem in final_list:
                if elem.equals(list1[i]):
                    has_equals = True
                    
            if not has_equals:
                final_list.append(list1[i])
                
        return final_list
    
    
    def maximize(self, student, focus, int_type):
        itvl_type = 0
        
        if int_type == self.high_share:
            itvl_type = 2
        elif int_type == self.low_share:
            itvl_type = 1
            
        intervals = self.manager.get_all_type(itvl_type)
        student_list = student.numbers
        subtrahend = 0
        focus_count = 0
        
        for i in range(len(student_list)):
            if student_list[i] != focus:
                subtrahend += intervals[student_list[i] - 1].high
            else:
                focus_count += 1
        
        for i in range(1, focus_count):
            subtrahend += intervals[focus - 1].high
            
        new_val = self.share.numerator - subtrahend
        
        if new_val < intervals[0].low:
            return -1
        else:
            return new_val
        
    
    def minimize(self, student, focus, int_type):
        itvl_type = 0
        
        if int_type == self.high_share:
            itvl_type = 2
        elif int_type == self.low_share:
            itvl_type = 1
            
        intervals = self.manager.get_all_type(itvl_type)
        student_list = student.numbers
        subtrahend = 0
        focus_count = 0
        
        for i in range(len(student_list)):
            if student_list[i] != focus:
                subtrahend += intervals[student_list[i] - 1].low
            else:
                focus_count += 1
            
        for i in range(1, focus_count):
            subtrahend += intervals[focus - 1].low
            
        new_val = intervals[focus - 1].low
        
        while (self.share.numerator - subtrahend - new_val) > 0:
            new_val += 0.5
            
        if new_val > intervals[len(intervals) - 1].high:
            return -1
        else:
            return new_val
        
    def refine(self, lst, int_type):
        itvl_type = 0
        if int_type == self.high_share:
            itvl_type = 2
        elif int_type == self.low_share:
            itvl_type = 1

        intervals = self.manager.get_all_type(itvl_type)
        to_be_added = []
        added = False

        for j in range(1, len(intervals) + 1):
            lowest_low = 0
            highest_high = intervals[j-1].high
            invalid = False

            for i in range(len(lst)):
                dmax = self.maximize(lst[i], j, int_type)
                dmin = self.minimize(lst[i], j, int_type)
                
                if lst[i].contains(j):
                    if self.max_or_min(lst[i], itvl_type, j) == 1:
                        if dmin != -1 and dmin <= highest_high and (dmin > lowest_low or lowest_low == 0):
                            lowest_low = dmin
                        elif dmin != -1 and dmin > highest_high and (dmax != -1 and (dmax < highest_high or highest_high == intervals[j-1].high)):
                            highest_high = dmax
                        elif dmin == -1:
                            invalid = True
                    else:
                        if dmax != -1 and dmax >= lowest_low and (dmax < highest_high or highest_high == intervals[j-1].high):
                            highest_high = dmax
                        elif dmax != -1 and dmax < lowest_low and (dmin != -1 and (dmin > lowest_low or lowest_low == 0)):
                            lowest_low = dmin
                        elif dmax == -1:
                            invalid = True

            if (lowest_low != intervals[j-1].low or highest_high != intervals[j-1].high) and intervals[j-1].contains_bound(lowest_low) and intervals[j-1].contains_bound(highest_high) and not invalid:
                if highest_high > lowest_low:
                    to_be_added.append(Interval(lowest_low, highest_high, 0))

        for elem in to_be_added:
            elem.low = self.myint(elem.low)
            elem.high = self.myint(elem.high)
            print("There are no shares in: (" + str(elem.low) + ", " + str(elem.high) + ")")           

            self.manager.add_interval(Interval(elem.low, elem.high, 0))
            self.buddy_one(Interval(elem.low, elem.high, 0))
            index_itvl = self.manager.find_interval(Interval(elem.low, elem.high, 0))
            index_itvl_buddy = self.manager.find_interval(Interval(self.share.denominator - elem.high, self.share.denominator - elem.low, 0))

            if self.manager.get(index_itvl + 1).type == 1 and self.low_share == 2 or self.manager.get(index_itvl - 1).type == 1 and self.low_share == 2:
                self.match_gap(self.manager.get(index_itvl))
            elif self.manager.get(index_itvl_buddy + 1).type == 1 and self.low_share == 2 or self.manager.get(index_itvl_buddy - 1).type == 1 and self.low_share == 2:
                self.match_gap(self.manager.get(index_itvl_buddy))
                
            added = True
            
        if not added:
            print("No gaps were found.")
            self.is_more_gaps = False

        self.easy_erik += 1
    
    def buddy_one(self, itvl):
        has_buddy = False
        for j in range(len(self.manager.get_all())):
            if self.manager.get(j).has_bound(self.share.denominator - itvl.low) and self.manager.get(j).has_bound(self.share.denominator - itvl.high):
                has_buddy = True

        if not has_buddy:
            bound1 = self.share.denominator - itvl.high
            bound2 = self.share.denominator - itvl.low
            
            print("There are no shares in (" + str(bound1)  + ", " + str(bound2) + ") as a result of buddying.")
   
            self.manager.add_interval(Interval(self.share.denominator - itvl.high, self.share.denominator - itvl.low, itvl.type, itvl.num_shares))
    
    def high_interval_students(self, num_high_interval):
        list_of_high_share = self.combination(num_high_interval, self.high_share)
        high_list = self.test(list_of_high_share, self.high_share)
        students = []
        temp = []

        i = 0
        while i <= len(high_list) - self.high_share:
            for j in range(i, i + self.high_share):
                temp.append(high_list[j])

            students.append(Student(temp, self.high_share))
            temp = []
            i += self.high_share

        return students
 
    def low_interval_students(self, num_low_interval):
        list_of_low_share = self.combination(num_low_interval, self.low_share)
        low_list = self.test(list_of_low_share, self.low_share)
        students = []
        temp = []

        i = 0
        while i <= len(low_list) - self.low_share:
            for j in range(i, i + self.low_share):
                temp.append(low_list[j])

            students.append(Student(temp, self.low_share))
            temp = []
            i += self.low_share

        return students
    
    def calc_num_high(self):
        num_high_interval = 0
        for elem in self.manager.get_all():
            if elem.type == 2:
                num_high_interval += 1

        return num_high_interval

    def calc_num_low(self):
        num_low_interval = 0
        for elem in self.manager.get_all():
            if elem.type == 1:
                num_low_interval += 1

        return num_low_interval

    def max_or_min(self, s, int_type, focus):
        new_type = self.high_share if int_type == 2 else self.low_share

        lst = self.manager.get_all_type(int_type)

        itvl = lst[focus - 1]

        student_list = s.numbers 

        dmax = 0
        dmin = 0

        for i in range(len(student_list)):
            dmax += lst[student_list[i] - 1].high 
            dmin += lst[student_list[i] - 1].low

        if abs(dmax - self.share.numerator) > abs(dmin - self.share.numerator) and itvl.contains_bound(self.minimize(s, focus, new_type)):
            return 1
        else:
            return 2
    
    def find_buddy_match_definitions(self, full_list, type_list, itvl):
        buddy_in_match = self.manager.get(self.find_buddy_definitions(full_list, itvl))
        match_of_buddy_index = self.manager.find_interval(Interval(self.share.numerator - buddy_in_match.high, self.share.numerator - buddy_in_match.low, buddy_in_match.type))

        if match_of_buddy_index != -1:
            match_of_buddy = self.manager.get(match_of_buddy_index)
            index = self.find_buddy_definitions(type_list, match_of_buddy)
            return index
        else:
            return -1
    
    def definitions(self, d_type):
        if d_type == 2:
            lst = self.manager.get_all_type(d_type)
            print("Definitions: All intervals are of " + str(self.high_share) + "-shares")
            union_list = []
            lone_share_list = []
            known_share_num = 0

            for i in range(1, len(lst) + 1):
                if lst[i - 1].num_shares == 0:
                    union_list.append(i)
                else:
                    lone_share_list.append(i) 
                    known_share_num += lst[i - 1].num_shares

            highest = union_list[len(union_list) - 1] if len(union_list) > 0 else 0
            highest_lone = lone_share_list[len(lone_share_list) - 1] if len(lone_share_list) > 0 else 0

            for i in range(1, len(lst) + 1):
                index = self.find_buddy_definitions(lst, lst[i - 1])
                index_match = self.find_buddy_match_definitions(self.manager.get_all(), lst, lst[i-1])

                if index != -1 and index + 1 < i:
                    print("J_" + str(i) + ": " + str(lst[i-1]), end =" ")
                    print("|J_" + str(index +1) + "| = |J_" + str(i) + "| ", end="")
                elif index_match != -1 and index_match + 1 < i:
                    print("J_" + str(i) + ": " + str(lst[i-1]), end=" ")
                    print("|J_" + str(index_match +1) + "| = |J_" + str(i) + "| ", end="")
                else:
                    print("J_" + str(i) + ": " + str(lst[i-1]), end=" ")

                if i == highest:
                    print("|",end="")
                    for j in range(len(union_list) - 1):
                        print("J_%d U " % union_list[j], end="")
                    f = union_list[len(union_list) - 1]
                    s = (self.num_low_share * self.low_share) + (self.num_high_share * self.high_share) - (2 * known_share_num)
                    print("J_{}| = {}".format(f, s))
                elif i == highest_lone:
                    print("|",end="")
                    for j in range(len(lone_share_list) - 1):
                        print("J_" + str(lone_share_list[j]) + " U ", end="")
                    print("J_" + str(lone_share_list[len(lone_share_list) - 1]) + "| = " + str(known_share_num))
                else:
                    print()
        else: 
            lst = self.manager.get_all_type(d_type)
            print("Definitions: All intervals are of " + str(self.low_share) + "-shares")
            union_list = []
            lone_share_list = []
            known_share_num = 0
            for i in range(1, len(lst) + 1):
                if lst[i - 1].num_shares == 0:
                    union_list.append(i)
                else:
                    lone_share_list.append(i)
                    known_share_num += lst[i - 1].num_shares

            highest = union_list[len(union_list) - 1] if len(union_list) > 0 else 0
            highest_lone = lone_share_list[len(lone_share_list) - 1] if len(lone_share_list) > 0 else 0

            for i in range(1, len(lst) + 1):
                index = self.find_buddy_definitions(lst, lst[i - 1])
                index_match = self.find_buddy_match_definitions(self.manager.get_all(), lst, lst[i - 1])

                if index != -1 and index + 1 < i:
                    print("I_" + str(i) + ": " + str(lst[i-1]), end=" ")
                    print("|I_" + str(index +1) + "| = |I_" + str(i) + "| ", end="")
                elif index_match != -1 and index_match + 1 < i:
                    print("I_" + str(i) + ": " + str(lst[i-1]), end=" ")
                else:
                    print("I_" + str(i) + ": " + str(lst[i-1]), end = " ")

                if i == highest:
                    print("|", end="")
                    for j in range(len(union_list) - 1):
                        print("I_" + str(union_list[j]) + " U ", end="")
                    print("I_" + str(union_list[len(union_list) - 1]) + "| = " + str((self.num_low_share * self.low_share) + (self.num_high_share * self.high_share) - (2 * known_share_num)))
                elif i == highest_lone:
                    print("|", end="")
                    for j in range(len(lone_share_list) - 1):
                        print("I_" + str(lone_share_list[j]) + " U ", end="")
                    print("I_" + str(lone_share_list[len(lone_share_list) - 1]) + "| = " + str(known_share_num))
                else:
                    print()

    
    def find_buddy_definitions(self, lst, itvl):
        index = -1
        for i in range(len(lst)):
            if lst[i].equals(Interval(self.share.denominator - itvl.high, self.share.denominator - itvl.low, itvl.type)):
                index = i

        return index

    
    def new_intervals(self):
        num_high_interval = self.calc_num_high()
        num_low_interval = self.calc_num_low()

        new_student_list = []
        if num_high_interval > 1 and num_low_interval == 1 or num_high_interval > num_low_interval:
            students = self.high_interval_students(num_high_interval)
            new_student_list = self.find_duplicates(students)   
            self.definitions(2)
            print("List of possible students: " + str(new_student_list))
            self.refine(new_student_list, self.high_share)
        elif num_low_interval > 1 and num_high_interval == 1 or num_low_interval > num_high_interval:
            students = self.low_interval_students(num_low_interval)
            new_student_list = self.find_duplicates(students)
            self.definitions(1)
            print("List of possible students: " + str(new_student_list))
            self.refine(new_student_list, self.low_share)

        return new_student_list
    
    def new_intervals2(self):
        num_high_interval = self.calc_num_high()
        num_low_interval = self.calc_num_low()

        new_student_list = []
        if num_high_interval > 1 and num_low_interval > num_high_interval:
            students = self.high_interval_students(num_high_interval)
            new_student_list = self.find_duplicates(students)
            self.definitions(2)
            print("List of possible students: " + str(new_student_list))
            self.refine(new_student_list, self.high_share)
        elif num_low_interval > 1 and num_high_interval > num_low_interval:
            students = self.low_interval_students(num_low_interval)
            new_student_list = self.find_duplicates(students)
            self.definitions(1)
            print("List of possible students: " + str(new_student_list))
            self.refine(new_student_list, self.low_share)

        return new_student_list
    
    def final_share_buddy(self):
        int_type = 0

        if self.calc_num_high() == 1:
            int_type = 2
        elif self.calc_num_low() == 1:
            int_type = 1
        elif self.calc_num_low() > self.calc_num_high():
            int_type = 2
        elif self.calc_num_high() > self.calc_num_low():
            int_type = 1

        if len(self.manager.get_all_type(int_type)) != 0:
            itvl = self.manager.get_all_type(int_type)[0]
            index = self.manager.find_buddy(itvl, self.share.denominator)
            if index != 1:
                self.manager.get(index).num_shares = itvl.num_shares
            else:
                print("BUDDYING ERROR")

        print()

    
    def find_initial_equations(self, int_type):
        interval_list = self.manager.get_all_type(int_type)
        final_list = []

        lone_share_list = []
        known_share_num = 0

        for i in range(1, len(interval_list) + 1):
            if interval_list[i - 1].num_shares != 0:
                lone_share_list.append(i)
                known_share_num += interval_list[i - 1].num_shares

        highest_lone = 0
        if len(lone_share_list) > 0:
            highest_lone = lone_share_list[len(lone_share_list) - 1]

        var = "J" if int_type == 2 else "I"

        for i in range(1, len(interval_list) + 1): 
            index = self.find_buddy_definitions(interval_list, interval_list[i - 1])
            index_match = self.find_buddy_match_definitions(self.manager.get_all(), interval_list, interval_list[i-1])

            if index != -1 and index + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index + 1))) 
                list2.append(Term(1, var + str(i)))

                final_list.append(Equation(list1, list2))
            elif index_match != -1 and index_match + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index_match + 1)))
                list2.append(Term(1, var + str(i)))

                final_list.append(Equation(list1, list2))

            if i == highest_lone:
                list1 = []
                list2 = []

                for j in range(len(lone_share_list)):
                    list1.append(Term(1, var + str(lone_share_list[j])))

                list2.append(Term(1, known_share_num))

                final_list.append(Equation(list1, list2))

        list1 = []
        list2 = []

        list1.append(Term(1, "s"))

        if int_type == 1:
            list2.append(Term(1, self.num_low_share))
        else:
            list2.append(Term(1, self.num_high_share))

        final_list.append(Equation(list1, list2))

        return final_list

    def find_initial_equations2(self, int_type):
        interval_list = self.manager.get_all_type(int_type)

        other_type = 1 if int_type == 2 else 2

        final_list = []
        lone_share_list = []

        known_share_num = 0

        for i in range(1, len(interval_list) + 1):
            if interval_list[i - 1].num_shares != 0:
                lone_share_list.append(i)
                known_share_num += interval_list[i - 1].num_shares

        highest_lone = 0

        if len(lone_share_list) > 0:
            highest_lone = lone_share_list[len(lone_share_list) - 1]

        var = "J" if int_type == 2 else "I"

        for i in range(1, len(interval_list) + 1):
            index = self.find_buddy_definitions(interval_list, interval_list[i - 1])
            index_match = self.find_buddy_match_definitions(self.manager.get_all(), interval_list, interval_list[i - 1])
            
            if index != -1 and index + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index + 1)))
                list2.append(Term(1, var + str(i)))

                final_list.append(Equation(list1, list2))
            elif index_match != -1 and index_match + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index_match + 1)))
                list2.append(Term(1, var + str(i)))

                final_list.append(Equation(list1, list2))
            if i == highest_lone:
                list1 = []
                list2 = []

                for j in range(len(lone_share_list)):
                    list1.append(Term(1, var + str(lone_share_list[j])))

                list2.append(Term(1, known_share_num))
                final_list.append(Equation(list1, list2))


        lone_share_list = []
        interval_list = self.manager.get_all_type(other_type)

        var = "J" if other_type == 2 else "I"

        known_share_num = 0

        for i in range(1, len(interval_list) + 1):
            if interval_list[i - 1].num_shares != 0:
                lone_share_list.append(i)
                known_share_num += interval_list[i - 1].num_shares

        highest_lone = 0

        if len(lone_share_list) > 0:
            highest_lone = lone_share_list[len(lone_share_list) - 1]

        for i in range(1, len(interval_list) + 1):
            index = self.find_buddy_definitions(interval_list, interval_list[i - 1])
            index_match = self.find_buddy_match_definitions(self.manager.get_all(), interval_list, interval_list[i - 1])
            
            if index != -1 and index + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index + 1)))
                list2.append(Term(1, var + str(i)))

                final_list.append(Equation(list1, list2))
            elif index_match != -1 and index_match + 1 < i:
                list1 = []
                list2 = []

                list1.append(Term(1, var + str(index_match + 1)))
                list2.append(Term(1, var + str(i)))

            if i == highest_lone:
                list1 = []
                list2 = []

                for j in range(len(lone_share_list)):
                    list1.append(Term(1, var + str(lone_share_list[j])))

                list2.append(Term(1, known_share_num))
                final_list.append(Equation(list1, list2))

        list1 = []
        list2 = []

        list1.append(Term(1, "sL"))
        list2.append(Term(1, self.num_low_share))

        final_list.append(Equation(list1, list2))

        list1 = []
        list2 = []

        list1.append(Term(1, "sH"))
        list2.append(Term(1, self.num_high_share))

        final_list.append(Equation(list1, list2))

        for i in range(len(self.manager.get_all_not(0))):
            index = self.find_buddy_definitions(self.manager.get_all_not(0), self.manager.get_all_not(0)[i])

            if index != -1 and self.find_buddy_definitions(self.manager.get_all_not(0), self.manager.get_all_not(0)[index]) < index and (not self.manager.get_all_not(0)[i] in self.manager.get_all_type(int_type) or not self.manager.get_all_not(0)[index] in self.manager.get_all_type(int_type)):
                list3 = []
                list4 = []

                var1 = ""
                var2 = ""

                true_index1 = 0
                true_index2 = 0

                if self.manager.get_all_not(0)[i].type == 1:
                    var1 += "I"
                    arr = self.manager.get_all_type(1)
                    true_index1 = arr.index(self.manager.get_all_not(0)[i])
                else:
                    var1 += "J"
                    arr = self.manager.get_all_type(2)
                    true_index1 = arr.index(self.manager.get_all_not(0)[i])

                if self.manager.get_all_not(0)[index].type == 1:
                    var2 += "I"
                    arr = self.manager.get_all_type(1)
                    true_index2 = arr.index(self.manager.get_all_not(0)[index])
                else:
                    var2 += "J"
                    arr = self.manager.get_all_type(2)
                    true_index2 = arr.index(self.manager.get_all_not(0)[index])

                list3.append(Term(1, var2 + str(true_index2 + 1)))
                list4.append(Term(1, var1 + str(true_index1 + 1)))
                final_list.append(Equation(list3, list4))

        return final_list
    
    def match_gap(self, itvl):
        numer = self.share.numerator 
        self.manager.add_interval(Interval(numer - itvl.high, numer - itvl.low, 0))
        bound1 = numer - itvl.high
        bound2 = numer - itvl.low
        
        print("There are no shares in (" + str(bound1)  + ", " + str(bound2) + ") as a result of match.")
        
        self.buddy_one(Interval(numer - itvl.high, numer - itvl.low, 0))
    
    def derive_equations(self, lst, initial_list):
        final_list = []

        for elem in initial_list:
            elem_arg1 = str(elem.arg1[0].variable)
            elem_arg2 = str(elem.arg2[0].variable)

            if elem_arg1[:1] == elem_arg2[:1] == "I":
                interval_num1 = int(elem_arg1[1:])
                arg1 = []

                for j in lst:
                    if j.contains(interval_num1) and len(j.numbers) == self.low_share:
                        arg1.append(Term(j.num_occurences(interval_num1), j))

                interval_num2 = int(elem_arg2[1:])
                arg2 = []

                for j in lst:
                    if j.contains(interval_num2) and len(j.numbers) == self.low_share:
                        arg2.append(Term(j.num_occurences(interval_num2), j))

                if len(arg2) == 0:
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))
            elif elem_arg1[:1] == elem_arg2[:1] == "J":
                interval_num1 = int(elem_arg1[1:])
                arg1 = []

                for j in lst:
                    if j.contains(interval_num1) and len(j.numbers) == self.high_share:
                        arg1.append(Term(j.num_occurences(interval_num1), j))

                interval_num2 = int(elem_arg2[1:])
                arg2 = []

                for j in lst:
                    if j.contains(interval_num2) and len(j.numbers) == self.high_share:
                        arg2.append(Term(j.num_occurences(interval_num2), j))

                if len(arg2) == 0:
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))
            elif elem_arg1[:1] == "I" and elem_arg2[:1] == "J":
                interval_num1 = int(elem_arg1[1:])
                arg1 = []

                for j in lst:
                    if j.contains(interval_num1) and len(j.numbers) == self.low_share:
                        arg1.append(Term(j.num_occurences(interval_num1), j))

                interval_num2 = int(elem_arg2[1:])
                arg2 = []

                for j in lst:
                    if j.contains(interval_num2) and len(j.numbers) == self.high_share:
                        arg2.append(Term(j.num_occurences(interval_num2), j))

                if len(arg2) == 0:
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))
            elif elem_arg1[:1] == "J" and elem_arg2[:1] == "I":
                interval_num1 = int(elem_arg1[1:])
                arg1 = []

                for j in lst:
                    if j.contains(interval_num1) and len(j.numbers) == self.high_share:
                        arg1.append(Term(j.num_occurences(interval_num1), j))

                interval_num2 = int(elem_arg2[1:])
                arg2 = []

                for j in lst:
                    if j.contains(interval_num2) and len(j.numbers) == self.low_share:
                        arg2.append(Term(j.num_occurences(interval_num2), j))

                if len(arg2) == 0:
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))
            elif elem_arg1[:1] == "I":
                arg1 = []

                for k in elem.arg1:
                    interval_num1 = int(k.variable[1:])

                    for j in lst:
                        if j.contains(interval_num1) and len(j.numbers) == self.low_share:
                            arg1.append(Term(j.num_occurences(interval_num1), j))

                arg2 = []
                arg2.append(Term(1, elem_arg2))

                final_list.append(Equation(arg1, arg2))
            elif elem_arg1[:1] == "J": 
                arg1 = []

                for k in elem.arg1:
                    interval_num1 = int(k.variable[1: len(k.variable)])

                    for j in lst:
                        if j.contains(interval_num1) and len(j.numbers) == self.high_share:
                            arg1.append(Term(j.num_occurences(interval_num1), j))

                arg2 = []
                arg2.append(Term(1, elem_arg2))

                if (len(arg2) == 0):
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))
            elif len(elem_arg1) == 2 and (elem_arg1[:2] == "sL" or elem_arg1[:2] == "sH"):
                arg1 = []
                arg2 = []

                if elem_arg1 == "sL":
                    for s in lst:
                        if len(s.numbers) == self.low_share:
                            arg1.append(Term(1, s))

                    arg2.append(Term(1, elem_arg2))
                    if len(arg2) == 0:
                        arg2.append(Term(1, 0.0))
                else:
                    for s in lst:
                        if len(s.numbers) == self.high_share:
                            arg1.append(Term(1, s))

                    arg2.append(Term(1, elem_arg2))

                final_list.append(Equation(arg1, arg2))
            else:
                arg1 = []
                arg2 = []

                for s in lst:
                    arg1.append(Term(1, s))

                arg2.append(Term(1, elem_arg2))
                if len(arg2) == 0:
                    arg2.append(Term(1, 0.0))

                final_list.append(Equation(arg1, arg2))

        return final_list 
    
    def convert_arg1(self, initial_list, lst):
        final_list = []

        size = len(lst) + 1

        for i in initial_list:
            eq = [0]*size
            eq[0] = 0

            for j in i.arg1:
                place = self.find_student(lst, j.variable) + 1

                if place != 0:
                    eq[place] = j.coeff 

            final_list.append(eq)

        return final_list
    
    def convert_arg2(self, initial_list):
        final_list = []

        for i in initial_list:
            final_list.append(i.arg2[0].variable)

        return final_list
    
    def find_student(self, lst, s):
        index = -1

        for i in range(len(lst)):
            if lst[i].equals(s):
                index = i

        return index 
    
    def find_student_in_equation(self, arg1, s):
        index = -1

        for i in range(len(arg1)):
            if arg1[i].variable.equals(s):
                index = i

        return index

    def shift_one_side(self, initial_list):
        for elem in initial_list:
            if (len(elem.arg1) == 0 and type(elem.arg2[0].variable) == Student) or type(elem.arg1[0].variable) == Student and type(elem.arg2[0].variable) == Student: 
                counter = 0
                while counter < len(elem.arg2):
                    index = self.find_student_in_equation(elem.arg1, elem.arg2[counter].variable)
                    if index == -1:
                        elem.arg2[counter].coeff = (-1 * elem.arg2[counter].coeff)
                        elem.arg1.append(elem.arg2[counter])
                        del elem.arg2[counter]
                    else:
                        coeff_of_right = elem.arg2[counter].coeff
                        elem.arg1[index].coeff = elem.arg1[index].coeff - coeff_of_right
                        del elem.arg2[counter]

                elem.arg2.append(Term(1, 0.0))

        return initial_list
    
    def linear_equations(self, lst, int_type, turn):
        initial = []

        if turn == 1:
            initial = self.find_initial_equations(int_type)
        else:
            initial = self.find_initial_equations2(int_type)

        derived = self.derive_equations(lst, initial)

        self.shift_one_side(derived) 

        for i in range(len(lst)):
            print("Let " + "var[" + str(i) + "] = " + str(lst[i]))
        print()

        arg1 = self.convert_arg1(derived, lst)
        arg2 = self.convert_arg2(derived)

        print("Left-hand arguments: ")

        for elem in arg1:
            for ele in elem:
                print(str(ele*1.0) + "  ", end="")
            print()

        print("Right-hand arguments: ")
        arg2 = list(map(lambda x: int(x)*1.0, arg2))
        print(arg2)
        print()

        # Linear Equations
        A = arg1
        b = arg2

        mod = LpProblem('prob')
        variables = LpVariable.dicts(name='x', indexs=range(len(A[0])), lowBound=0, cat="Integer")
        for row, rhs in zip(A, b):
            mod += sum([row[i]*variables[i] for i in range(len(row))]) == rhs 
        
        mod.solve()
        print(LpStatus[mod.status])

        if not mod.objective:
            print("Value of objective function:", mod.objective)
        [print("Value of var["+ str(i) + "] = " + str(variables[i + 1].value())) for i in range(len(A[0]) - 1)]  
        var = [variables[i + 1].value() for i in range(len(A[0]) - 1)]

        x = [variables[i+1].value() for i in range(len(A[0]) - 1)]
        for arr in A:
            del arr[0]

        if (var == [0.0]*(len(A[0]) - 1) or self.contains_floats(var)):
            print("Solutions: 0")
            return 0
        elif LpStatus[mod.status] == "Infeasible":
            print("Solutions: 0")
            return 0
        else:
            print("Solutions: 1") 
            return 1

    def contains_floats(self, arr):
        flag = False
        for i in arr:
            if int(i) != i:
                flag = True

        return flag

    def buddy_two(self, itvl):
        has_buddy = False
        
        for j in range(len(self.manager.get_all())):
            if self.manager.get(j).has_bound(self.share.denominator - itvl.low) and self.manager.get(j).has_bound(self.share.denominator - itvl.high):
                has_buddy = True
                
        if not has_buddy:
            self.manager.add_interval(Interval(self.share.denominator - itvl.high, self.share.denominator - itvl.low, itvl.type, itvl.num_shares))
                     
    def match(self):
        midpoint = 0
        
        for i in range(len(self.manager.get_all())):
            new_num_shares = self.manager.get(i).num_shares
            
            if self.low_share == 2 and self.manager.get(i).type == 1 or self.high_share == 2 and self.manager.get(i).type == 2:
                if self.share.numerator - self.manager.get(i).low == self.manager.get(i).high:
                    low = self.manager.get(i).low
                    midpoint = self.myint((self.manager.get(i).high + self.manager.get(i).low)/2)
                    
                    self.manager.add_interval(Interval(low, midpoint, 2, self.myint(new_num_shares/2)))
                    index = self.manager.find_interval(Interval(low, midpoint, 2, self.myint(new_num_shares/2)))
                    self.manager.get(index + 1).num_shares = self.myint(new_num_shares/2)
                    
                    self.buddy_two(Interval(low, midpoint, 2, self.myint(new_num_shares/2))) 
                    index = self.manager.find_interval(Interval(self.share.denominator - midpoint, self.share.denominator - low, 2, self.myint(new_num_shares/2)))
                    self.manager.get(index - 1).num_shares = self.myint(new_num_shares/2)
                    
    def is_bill(self):
        isB = False
        midpoint = self.myint((self.manager.get_highest_bound() + self.manager.get_lowest_bound())/2)
        
        size = len(self.manager.get_all())
        if size == 2:
            if self.manager.get(size - 1).low >= midpoint:
                if self.manager.get(size - 1).type == 2 and self.high_share * self.num_high_share > self.m:
                    isB = True
                    print("All " + str(self.high_share) + " shares > than the midpoint (" + str(midpoint) + ")")
                    print("Contradiction. " + str(self.high_share) + " * s" + str(self.high_share) + " > m since " + str(self.high_share) + " * " + str(self.num_high_share) + " > " + str(self.m))
                elif self.manager.get(size - 1).type == 1 and self.low_share * self.num_low_share > self.m:
                    isB = True
                    print("All " + str(self.low_share) + " shares > than the midpoint (" + str(midpoint) + ")")
                    print("Contradiction. " + str(self.low_share) + " * s" + str(self.low_share) + " > m since " + str(self.low_share) + " * " + str(self.num_low_share) + " > " + str(self.m))
            elif self.manager.get(0).high <= midpoint:
                if self.manager.get(0).type == 2 and self.high_share * self.num_high_share > self.m:
                    isB = True
                    print("All " + str(self.high_share) + " shares < than the midpoint (" + str(midpoint) + ")")
                    print("Contradiction. " + str(self.high_share) + " * s" + str(self.high_share) + " > m since " + str(self.high_share) + " * " + str(self.num_high_share) + " > " + str(self.m))
                elif self.manager.get(0).type == 1 and self.low_share * self.num_low_share > self.m:
                    isB = True
                    print("All " + str(self.low_share) + " shares < than the midpoint (" + str(midpoint) + ")")
                    print("Contradiction. " + str(self.low_share) + " * s" + str(self.low_share) + " > m since " + str(self.low_share) + " * " + str(self.num_low_share) + " > " + str(self.m))
                    
        return isB    
          
    def main(self):
        try:
            self.three_share_test()
            self.set_interval()
            
            itvl = self.find_upper_range() 
            new_bound = self.find_lower_range()
            
            if len(self.manager.get_all()) == 2 and self.manager.get(1).type == 0 and itvl.low <= self.manager.get(1).low and itvl.high >= self.manager.get(1).high:
                if itvl.low < self.manager.get(0).low:
                    self.manager.remove(1)
                    self.manager.add_ignore_everything(0, itvl)
                elif self.manager.get(0).low < itvl.low:
                    self.manager.remove(1)
                    self.manager.add_ignore_everything(1, itvl)
                elif self.manager.get(1).high == itvl.high:
                    self.manager.remove(1)
                    self.manager.add_ignore_everything(1, itvl)
                else:
                    print("Whattt???")
            elif len(self.manager.get_all()) == 1:
                self.manager.add_ignore_everything(1, itvl)
            else:
                self.manager.add_ignore_share(itvl) 
                
            if len(self.manager.get_all_type(1)) > 0: 
                itvl_low = self.manager.get_all_type(1)[0]
                print("%s-share interval: %s" % (self.low_share, itvl_low))
            else:
                print("There are no intervals with %s-shares" % (self.low_share))
                
            if len(self.manager.get_all_type(2)) > 0:
                itvl_high = self.manager.get_all_type(2)[0]
                print("%s-share interval: %s" % (self.high_share, itvl_high))
            else:
                print("There are no intervals with %s-shares" % (self.high_share))

            self.find_blocks()
            try:
                if not self.is_bill(): 
                    if (len(self.manager.get_all()) == 2 and self.manager.get(0).high > self.manager.get(1).low):
                        print("Tried f(" + str(self.m) + "," + str(self.s) + "\\le " + str(self.alpha) + " and failed BOO-overlapping intervals") 
                    else:
                        self.mid_point()
                        self.final_share_buddy()
                        self.match() 
                        final_list = []
                        
                        while self.is_more_gaps:
                            final_list = self.new_intervals()
                            self.match()
                            
                        int_type = 0
                        
                        if self.calc_num_high() == 1:
                            int_type = 1
                        elif self.calc_num_low() == 1:
                            int_type = 2
                        elif self.calc_num_low() > self.calc_num_high():
                            int_type = 1
                        elif self.calc_num_high() > self.calc_num_low():
                            int_type = 2    
                        else:
                            print("EQUAL NUMBER OF INTERVALS OF BOTH SHARE NUMBERS?")
                        
                        num_solutions = -1
                        
                        if len(final_list) != 0:
                            num_solutions = self.linear_equations(final_list, int_type, self.turn)

                        if num_solutions == 1: 
                            self.is_more_gaps = True
                            num_high_interval = self.calc_num_high()
                            num_low_interval = self.calc_num_low()

                            if num_high_interval > 1 and num_low_interval > num_high_interval or num_low_interval > 1 and num_high_interval > num_low_interval:
                                new_list = []

                                while self.is_more_gaps:
                                    new_list = self.new_intervals2()

                                for elem in new_list:
                                    final_list.append(elem)

                                self.turn += 1
                                num_sol = self.linear_equations(final_list, int_type, self.turn)

                                if num_sol != 1:
                                    print()
                                    print("Tried f(%s, %s) \\le %s and succeeded YEAH-BILL2" % (self.m, self.s, self.alpha))
                            else:
                                print("Tried f(%s, %s) \\le %s and failed BOO-ERIK FAILED" % (self.m, self.s, self.alpha))
                        elif num_solutions == 0:
                            print()
                            if self.easy_erik == 1:
                                print("Tried f(%s, %s) \\le %s and succeeded YEAH-EERIK" % (self.m, self.s, self.alpha))
                            else:
                                print("Tried f(%s, %s) \\le %s and succeeded YEAH-ERIK" % (self.m, self.s, self.alpha))
                        elif num_solutions == -1:
                            self.is_more_gaps = True
                            num_high_interval = self.calc_num_high()
                            num_low_interval = self.calc_num_low()

                            if num_high_interval > 1 and num_low_interval > num_high_interval or num_low_interval > 1 and num_high_interval > num_low_interval:
                                new_list = []

                                while self.is_more_gaps:
                                    new_list = self.new_intervals2()

                                for elem in new_list:
                                    final_list.append(elem)

                                self.turn += 1

                                num_sol = self.linear_equations(final_list, int_type, self.turn)

                                if num_sol != 1:
                                    print("WARNING: Solved by bill2 though ERIK did not work for the first set of shares")

                            else:
                                print("Bill 2 was attempted because there were 0 students that used the first number-of-share and subsequently failed.")
                                print()
                else:
                    raise BillCaseException("BILL-CASE")

            except BillCaseException:
                print("Tried f(%s, %s)\\le %s and succeeded YEAH-BILL1" % (self.m, self.s, self.alpha))  
            
        except ThreeShareException:
            print("Tried f(%s, %s) \\le %s and failed BOO-three numbers of shares" % (self.m, self.s, self.alpha))  
            print()
        
default_stdout = sys.stdout
file_handle = open('output.txt', 'w')
try:
    sys.stdout = file_handle 
    args = sys.argv

    if len(args) == 5:
        muffin = MuffinsAuto(int(args[1]), int(args[2]), int(args[3]), int(args[4]))
        muffin.main()
    elif len(args) == 2:
        with open(args[1], "r") as ins:
            arr = []
            for line in ins:
                arr.append(line.rstrip('\n').split())

        for elem in arr:
            name = elem[0]+"-"+elem[1]+'.txt'
            file = open(name, 'w')
            sys.stdout = file 

            muffin = MuffinsAuto(int(elem[0]), int(elem[1]), int(elem[2]), int(elem[3]))
            muffin.main()

            file.close()

            with open(name, 'r') as output:
                data = output.read()

            if "EERIK" in data:
                os.rename(name, elem[0]+"-"+elem[1]+'.EERIK.txt')
            elif "fail" in data:
                os.rename(name, elem[0]+"-"+elem[1]+'.OPEN.txt')
            elif "ERIK" in data:
                os.rename(name, elem[0]+"-"+elem[1]+'.ERIK.txt')
            elif "BILL1" in data:
                os.rename(name, elem[0]+"-"+elem[1]+'.BILL1.txt')
            elif "BILL2" in data:
                os.rename(name, elem[0]+"-"+elem[1]+'.BILL2.txt')
finally:
    # Restores stdout, even if an exception occurs.
    sys.stdout = default_stdout
