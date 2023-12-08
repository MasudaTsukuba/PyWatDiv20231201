# volatility_gen.py
# Python version of WatDiv/volatility_gen
# T. Masuda, 2023/12/1
import bisect
import random
import sys
from enum import Enum, auto
import math
from scipy.special import erfinv

# #include "volatility_gen.h"

# #include <boost/lexical_cast.hpp>
# #include <boost/math/constants/constants.hpp>
# #include <boost/math/special_functions/erf.hpp>

# #include <algorithm>
# #include <cmath>
# #include <fstream>
# #include <iostream>

# using namespace std;

# namespace VOLATILITY_MODEL {
#     enum enum_t {
#         NORMAL_DIST,
#         LAPLACE_DIST,
#         LOGISTIC_DIST,
#         CAUCHY_DIST,
#         UNDEFINED
#     };
# };


class VolatilityModel(Enum):
    NORMAL_DIST = auto
    LAPLACE_DIST = auto
    LOGISTIC_DIST = auto
    CAUCHY_DIST = auto
    UNDEFINED = auto


def sgn(val):
    """

    Args:

    Returns:

    """
    return (0 < val) - (val < 0)


class VolatilityGen:
    """

    Attributes:

    """

    def __init__(self, frequency_sample_file=None, location=None, distribution=None, volatility_sample_file=None, rhs=None):
        """

        Args:

        Returns:

        """
        #         mt19937 * _rd_gen;
        self._rd_gen = None
        #         uniform_real_distribution<> * _uniform_real_distribution;
        self._uniform_real_distribution = None

        #         bool _instantiated;
        self._instantiated: bool = False
        #         int _model_size;
        self._model_size: int = 0
        #         double _location;
        self._location: float = 0.0
        #         double _max_frequency;
        self._max_frequency: float = 0.0
        #         VOLATILITY_MODEL::enum_t _distribution;
        self._distribution = None
        #         vector<double> _frequency_sample;
        self._frequency_sample = []
        #         vector<double> _gen_frequency;
        self._gen_frequency = []
        #         vector<double> _volatility_sample;
        self._volatility_sample = []
        #         vector<double> _gen_volatility;
        self._gen_volatility = []

        #         double * _dynamic_frequency_array;
        self._dynamic_frequency_array = []
        #         double * _cum_density_array;
        self._cum_density_array = []

        #         template <typename T> int sgn(T val) const {
        #             return (T(0) < val) - (val < T(0));
        #         }
        pass

# volatility_gen::volatility_gen(const char * frequency_sample_file, double location, VOLATILITY_MODEL::enum_t distribution, const char * volatility_sample_file){
        if frequency_sample_file is not None and location is not None and distribution is not None and volatility_sample_file is not None:
#     random_device rd;

#     _rd_gen = new mt19937 (rd());

#     _uniform_real_distribution = new uniform_real_distribution<> (0, 1);


#     _instantiated = false;
            self._instantiated = False
#     _model_size = 0;
            self._model_size = 0

#     _location = location;
            self._location = location
#     _distribution = distribution;
            self._distribution = distribution

#     string token;
#     ifstream ifs_freq (frequency_sample_file);
            with open(frequency_sample_file, 'r') as ifs_freq:
                lines = ifs_freq.readlines()
#     if (!ifs_freq){
#         cerr << "[volatility_gen::volatility_gen]\tFile "<<frequency_sample_file<<" does not exist..." << "\n";
#         exit(0);
#     }
#     while (ifs_freq>>token){
                for token in lines:
#         double freq_value = boost::lexical_cast<double>(token);
                    freq_value = float(token)
#         _frequency_sample.push_back(freq_value);
                    self._frequency_sample.append(freq_value)
#     }
#     ifs_freq.close();
#     sort(_frequency_sample.begin(), _frequency_sample.end());
            self._frequency_sample.sort()

#     bool zero_included = false;
            zero_included = False
#     ifstream ifs_vol (volatility_sample_file);
            with open(volatility_sample_file, 'r') as ifs_vol:
                lines = ifs_vol.readlines()
#     if (!ifs_vol){
#         cerr << "[volatility_gen::volatility_gen]\tFile "<<volatility_sample_file<<" does not exist..." << "\n";
#         exit(0);
#     }
#     while (ifs_vol>>token){
                for token in lines:
#         double vol_value = boost::lexical_cast<double>(token);
                    vol_value = float(token)
#         if (vol_value==0.0){
                    if vol_value == 0.0:
#             zero_included = true;
                        zero_included = True
#         }
#         _volatility_sample.push_back(vol_value);
                    self._volatility_sample.append(vol_value)
#     }
#     ifs_vol.close();
#     if (!zero_included){
                if not zero_included:
#         _volatility_sample.push_back(0.0);
                    self._volatility_sample.append(0.0)
#     }
#     sort(_volatility_sample.begin(), _volatility_sample.end());
            self._volatility_sample.sort()
# }
        pass  # end of volatility_gen::volatility_gen

# volatility_gen::volatility_gen(const volatility_gen & rhs){
        if rhs is not None:
#     random_device rd;
#     _rd_gen = new mt19937 (rd());
#     _uniform_real_distribution = new uniform_real_distribution<> (0, 1);

#     _instantiated = rhs._instantiated;
            self._instantiated = rhs._instantiated
#     _model_size = rhs._model_size;
            self._model_size = rhs._model_size
#     _location = rhs._location;
            self._location = rhs._location
#     _max_frequency = rhs._max_frequency;
            self._max_frequency = rhs._max_frequency
#     _distribution = rhs._distribution;
            self._distribution = rhs._distribution
#     _frequency_sample.insert(_frequency_sample.begin(), rhs._frequency_sample.cbegin(), rhs._frequency_sample.cend());
            self._frequency_sample += rhs._frequency_sample
#     _gen_frequency.insert(_gen_frequency.begin(), rhs._gen_frequency.cbegin(), rhs._gen_frequency.cend());
            self._gen_frequency += rhs._gen_frequency
#     _volatility_sample.insert(_volatility_sample.begin(), rhs._volatility_sample.cbegin(), rhs._volatility_sample.cend());
            self._volatility_sample += rhs._volatility_sample
#     _gen_volatility.insert(_gen_volatility.begin(), rhs._gen_volatility.cbegin(), rhs._gen_volatility.cend());
            self._gen_volatility += rhs._gen_volatility

#     _dynamic_frequency_array = new double [_model_size];
            self._dynamic_frequency_array = []
#     _cum_density_array = new double [_model_size];
            self._cum_density_array = []
#     for (int i=0; i<_model_size; i++){
            for i in range(self._model_size):
#         _dynamic_frequency_array[i] = rhs._dynamic_frequency_array[i];
                self._dynamic_frequency_array.append(rhs._dynamic_frequency_array[i])
#         _cum_density_array[i] = rhs._cum_density_array[i];
                self._cum_density_array.append(rhs._cum_density_array[i])
#     }
# }
        pass  # end of volatility_gen::volatility_gen(const volatility_gen & rhs)

# volatility_gen::~volatility_gen(){
#     delete _rd_gen;
#     delete _uniform_real_distribution;
#     if (_model_size>0){
#         delete [] _cum_density_array;
#         delete [] _dynamic_frequency_array;
#     }
# }

# void volatility_gen::initialize (int model_size){
    def initialize(self, model_size):
        """

        Args:

        Returns:

        """
#     _model_size = model_size;
        self._model_size = model_size
#     _dynamic_frequency_array = new double [_model_size];
        self._dynamic_frequency_array = []
#     _cum_density_array = new double[_model_size];
        self._cum_density_array = []
#     _max_frequency = 0.0;
        self._max_frequency = 0.0

#     double total_frequency = 0.0;
        total_frequency: float = 0.0
#     for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         double rand_index = next_uniform() * (double)(_frequency_sample.size()-1);
            rand_index = self.next_uniform() * (len(self._frequency_sample) - 1)

#         double range_min = _frequency_sample[floor(rand_index)];
            range_min = self._frequency_sample[math.floor(rand_index)]
#         double range_max = _frequency_sample[ceil(rand_index)];
            range_max = self._frequency_sample[math.ceil(rand_index)]
#         double interpolated = range_min + ((range_max - range_min) * (rand_index - floor(rand_index)));
            interpolated = range_min + ((range_max - range_min) * (rand_index - math.floor(rand_index)))
#         if (interpolated > _max_frequency){
            if interpolated > self._max_frequency:
#             _max_frequency = interpolated;
                self._max_frequency = interpolated
#         }
#         _gen_frequency.push_back(interpolated);
            self._gen_frequency.append(interpolated)
#         _dynamic_frequency_array[i] = interpolated;
            self._dynamic_frequency_array[i] = interpolated
#         total_frequency = total_frequency + interpolated;
            total_frequency = total_frequency + interpolated
#     }
            pass  # end of for
#     for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         double rand_index = next_uniform() * (double)(_volatility_sample.size()-1);
            rand_index = self.next_uniform() * (len(self._volatility_sample) - 1)
#         double range_min = _volatility_sample[floor(rand_index)];
            range_min = self._volatility_sample[math.floor(rand_index)]
#         double range_max = _volatility_sample[ceil(rand_index)];
            range_max = self._volatility_sample[math.ceil(rand_index)]
#         double interpolated = range_min + ((range_max - range_min) * (rand_index - floor(rand_index)));
            interpolated = range_min + ((range_max - range_min) * (rand_index - math.floor(rand_index)))
#         _gen_volatility.push_back(interpolated);
            self._gen_volatility.append(interpolated)
#     }
            pass  # end of for

#     for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         if (i==0){
            if i == 0:
#             _cum_density_array[i] = _dynamic_frequency_array[i] / total_frequency;
                self._cum_density_array[i] = self._dynamic_frequency_array[i] / total_frequency
#         } else {
            else:
#             _cum_density_array[i] = _cum_density_array[i-1] + (_dynamic_frequency_array[i] / total_frequency);
                self._cum_density_array[i] = self._cum_density_array[i - 1] + (self._dynamic_frequency_array[i] / total_frequency)
#         }
#     }
            pass  # end of for
#     _max_frequency = _max_frequency * 1000;
        self._max_frequency = self._max_frequency * 1000
#     _instantiated = true;
        self._instantiated = True
#
#     /// Some debugging...
#     //for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         //cout << "Volatility " << i << "=" << _gen_volatility[i] << "\n";
            print(f"Volatility {i} ={self._gen_volatility[i]}")
#     //}
            pass  # end of for
#     warmup();
        self.warmup()
# }
        pass  # end of volatility_gen::initialize

# bool volatility_gen::is_initialized () const{
    def is_initialized(self):
        """

        Args:

        Returns:

        """
#     return _instantiated;
        return self._instantiated
# }
        pass  # end of volatility_gen::is_initialized

# double volatility_gen::advance (){
    def advance(self):
        """

        Args:

        Returns:

        """
#     double total_frequency = 0.0;
        total_frequency = 0.0
#     for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         double ln_inc = get_next_ln_increase(i);
            ln_inc = self.get_next_ln_increase(i)
#         _dynamic_frequency_array[i] = _dynamic_frequency_array[i] * exp(ln_inc);
            self._dynamic_frequency_array[i] = self._dynamic_frequency_array[i] * math.exp(ln_inc)
#         if (_dynamic_frequency_array[i] < 1.0){
            if self._dynamic_frequency_array[i] < 1.0:
#             _dynamic_frequency_array[i] = 1.0;
                self._dynamic_frequency_array[i] = 1.0
#         }
#         if (_dynamic_frequency_array[i] > _max_frequency){
            if self._dynamic_frequency_array[i] > self._max_frequency:
#             _dynamic_frequency_array[i] = _max_frequency;
                self._dynamic_frequency_array[i] = self._max_frequency
#         }
#         total_frequency += _dynamic_frequency_array[i];
            total_frequency += self._dynamic_frequency_array[i]
#         //cout << "log-increase [" << i << "] = " << ln_inc << "\n";
#         //cout << "ovr-increase [" << i << "] = " << exp(ln_inc) << "\n";
#         //cout << "frequency [" << i << "] = " << _dynamic_frequency_array[i] << "\n";
#     }
            pass  # end of for
#     for (int i=0; i<_model_size; i++){
        for i in range(self._model_size):
#         if (i==0){
            if i == 0:
#             _cum_density_array[i] = _dynamic_frequency_array[i] / total_frequency;
                self._cum_density_array[i] = self._dynamic_frequency_array[i] / total_frequency
#         } else {
            else:
#             _cum_density_array[i] = _cum_density_array[i-1] + (_dynamic_frequency_array[i] / total_frequency);
                self._cum_density_array[i] = self._cum_density_array[i - 1] + (self._dynamic_frequency_array[i] / total_frequency)
#         }
#         //cout << "probability [" << i << "] = " << (_dynamic_frequency_array[i] / total_frequency) << "\n";
#     }
            pass  # end of for
#     //cout << "Total frequency = " << total_frequency << "\n";
#     cerr << "Advancing..." << "\n";
#     return total_frequency;
        return total_frequency
# }
        pass  # end of volatility_gen::advance

# void volatility_gen::warmup (){
    def warmup(self):
        """

        Args:

        Returns:

        """
#     cerr << "[volatility_gen::warmup()]\tWarmup phase has started..." << "\n";

#     int iteration_limit = 1000;
        iteration_limit = 1000
#     int window_size = 5;
        window_size = 5
#     double termination_threshold = 0.01;
        termination_threshold = 0.01

#     vector<double> history;
        history = []
#     double prev_total = 1;
        prev_total = 1
#     for (int i=0; i<iteration_limit; i++){
        for i in range(iteration_limit):
#         double cur_total = advance();
            cur_total = self.advance()
#         double log_increase = log(cur_total / prev_total);
            log_increase = math.log(cur_total / prev_total)
#         prev_total = cur_total;
            prev_total = cur_total
#         history.push_back(log_increase);
            history.append(log_increase)
#         if (history.size()>=window_size){
            if len(history) >= window_size:
#             double sum = 0.0;
                sum0 = 0.0
#             for (int j=(history.size()-window_size); j<history.size(); j++){
                for j in range(len(history) - window_size, len(history)):
#                 sum += history[j];
                    sum0 += history[j]
#             }
                    pass  # end of for
#             double average = sum / ((double) window_size);
                average = sum0 / float(window_size)
#             if (average < termination_threshold){
                if average < termination_threshold:
#                 cerr << "[volatility_gen::warmup()]\tWarmup phase ended after " << i << " iterations..." << "\n";
                    print(f"[volatility_gen::warmup()]\tWarmup phase ended after {i} iterations...")
#                 break;
                    break
#             }
                    pass  # end of if
#         }
                pass  # end of if
#     }
            pass  # end of for
# }
        pass  # end of volatility_gen::warmup

# int volatility_gen::get_model_size () const{
    def get_model_size(self):
        """

        Args:

        Returns:

        """
#     return _gen_volatility.size();
        return len(self._gen_volatility)
# }
        pass  # end of volatility_gen::get_model_size

# double volatility_gen::get_initial_frequency_value (int index) const{
    def get_initial_frequency_value(self, index):
        """

        Args:

        Returns:

        """
#     return _gen_frequency[index];
        return self._gen_frequency[index]
# }
        pass  # end of volatility_gen::get_initial_frequency_value

# double volatility_gen::get_volatility_value (int index) const{
    def get_volatility_value(self, index):
        """

        Args:

        Returns:

        """
#     return _gen_volatility[index];
        return self._gen_frequency[index]
# }
        pass  # end of volatility_gen::get_volatility_value

# double volatility_gen::get_next_ln_increase (int index) const{
    def get_next_ln_increase(self, index):
        """

        Args:

        Returns:

        """
        #     double pr = next_uniform();
        pr: float = self.next_uniform()
#     double result = 0.0;
        result = 0.0
#     switch (_distribution){
#         case VOLATILITY_MODEL::NORMAL_DIST:{
        if self._distribution == VolatilityModel.NORMAL_DIST:
#             result = _location + (_gen_volatility[index] * sqrt(2.0) * boost::math::erf_inv((2 * pr) - 1.0));
            result = self._location + (self._gen_volatility[index] * math.sqrt(2.0) * erfinv((2.0 * pr) - 1.0))
#             break;
#         }
#         case VOLATILITY_MODEL::CAUCHY_DIST:{
        elif self._distribution == VolatilityModel.CAUCHY_DIST:
#             result = _location + (_gen_volatility[index] * tan(boost::math::constants::pi<double>() * (pr - 0.5)));
            result = self._location + (self._gen_volatility[index] * math.tan(math.pi * (pr - 0.5)))
#             break;
#         }
#         case VOLATILITY_MODEL::LAPLACE_DIST:{
        elif self._distribution == VolatilityModel.LAPLACE_DIST:
#             result = _location - (_gen_volatility[index] * sgn<double>(pr - 0.5) * log( 1 - (2 * abs(pr - 0.5)) ));
            result = self._location + (self._gen_volatility[index] * sgn(pr - 0.5) * math.log(1.0 - (2.0 * abs(pr - 0.5))))
#             break;
#         }
#         case VOLATILITY_MODEL::LOGISTIC_DIST:{
        elif self._distribution == VolatilityModel.LOGISTIC_DIST:
#             result = _location + (_gen_volatility[index] * log(pr/(1-pr)));
            result = self._location + (self._gen_volatility[index] * math.log(pr/(1.0 - pr)))
#             break;
#         }
#     }
        pass  # end of if
#     return result;
        return result
# }
        pass  # end of volatility_gen::get_next_ln_increase

# double volatility_gen::get_probability (int index) const{
    def get_probability(self, index):
        """

        Args:

        Returns:

        """
        #     if (index==0){
        if index == 0:
#         return _cum_density_array[index];
            return self._cum_density_array[index]
#     } else {
        else:
#         return (_cum_density_array[index] - _cum_density_array[(index-1)]);
            return self._cum_density_array[index] - self._cum_density_array[index - 1]
#     }
        pass  # end of if
# }
        pass  # end of volatility_gen::get_probability

# int volatility_gen::next_rand_index () const{
    def next_rand_index(self):
        """

        Args:

        Returns:

        """
        #     double pivot = next_uniform();
        pivot = self.next_uniform()
#     double * location = lower_bound(_cum_density_array, _cum_density_array+_model_size, pivot);
        location = bisect.bisect_left(self._cum_density_array, pivot)
#     return (location-_cum_density_array);
        return location - self._cum_density_array[0]
# }
        pass  # end of volatility_gen::next_rand_index

# double volatility_gen::next_uniform() const{
    def next_uniform(self):
        """

        Args:

        Returns:

        """
        #     return (*_uniform_real_distribution)(*_rd_gen);
        return random.uniform(0, 1)
# }
        pass  # end of volatility_gen::next_uniform

    @staticmethod
# volatility_gen * volatility_gen::parse (const string & line, string & name, float & advance_pr){
    def parse(line, name=None, advance_pr=None):
        """

        Args:

        Returns:

        """
        #     double location = 0.0;
        location = 0.0
#     VOLATILITY_MODEL::enum_t distribution = VOLATILITY_MODEL::UNDEFINED;
        distribution = VolatilityModel.UNDEFINED
#     string frequency_sample_file;
        frequency_sample_file = ""
#     string volatility_sample_file;
        volatility_sample_file = ""

#     stringstream parser(line);
        parser = line.split(' ')
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("#dynamic")!=0){
                if token != "#dynamic":
#                     cerr<<"[volatility_gen::parse()]\tExpecting #dynamic..."<<"\n";
                    print("[volatility_gen::parse()]\tExpecting #dynamic...")
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            elif index == 1:
#                 name = token;
                name = token
#                 break;
#             }
#             case 2:{
            elif index == 2:
#                 location = boost::lexical_cast<double>(token);
                location = float(token)
#                 break;
#             }
#             case 3:{
            elif index == 3:
#                 if (token.compare("NORMAL")==0 || token.compare("normal")==0){
                if token == "NORMAL" or token == "normal":
#                     distribution = VOLATILITY_MODEL::NORMAL_DIST;
                    distribution = VolatilityModel.NORMAL_DIST
#                 } else if (token.compare("CAUCHY")==0 || token.compare("cauchy")==0){
                elif token == "CAUCHY" or token == "cauchy":
#                     distribution = VOLATILITY_MODEL::CAUCHY_DIST;
                    distribution = VolatilityModel.CAUCHY_DIST
#                 } else if (token.compare("LAPLACE")==0 || token.compare("laplace")==0){
                elif token == "LAPLACE" or token == "laplace":
#                     distribution = VOLATILITY_MODEL::LAPLACE_DIST;
                    distribution = VolatilityModel.LAPLACE_DIST
#                 } else if (token.compare("LOGISTIC")==0 || token.compare("logistic")==0){
                elif token == "LOGISTIC" or token == "logistic":
#                     distribution = VOLATILITY_MODEL::LOGISTIC_DIST;
                    distribution = VolatilityModel.LOGISTIC_DIST
#                 } else {
                else:
#                     cerr<<"[volatility_gen::parse()]\tUnsupported distribution..."<<"\n";
                    print("[volatility_gen::parse()]\tUnsupported distribution...")
#                     cerr<<"                         \tExpecting one of NORMAL, CAUCHY, LAPLACE, LOGISTIC..."<<"\n";
                    print("                         \tExpecting one of NORMAL, CAUCHY, LAPLACE, LOGISTIC...")
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 4:{
            elif index == 4:
#                 frequency_sample_file = token;
                frequency_sample_file = token
#                 break;
#             }
#             case 5:{
            elif index == 5:
#                 volatility_sample_file = token;
                volatility_sample_file = token
#                 break;
#             }
#             case 6:{
            elif index == 6:
#                 advance_pr = boost::lexical_cast<float>(token);
                advance_pr = float(token)
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#     if (index!=7){
        if index != 7:
#         cerr<<"[volatility_gen::parse()]\tUnsupported number of arguments..."<<"\n";
            print("[volatility_gen::parse()]\tUnsupported number of arguments...")
#         exit(0);
            sys.exit(0)
#     }
#     return new volatility_gen(frequency_sample_file.c_str(), location, distribution, volatility_sample_file.c_str());
        return VolatilityGen(frequency_sample_file, location, distribution, volatility_sample_file)
# }
        pass  # end of volatility_gen::parse

    @staticmethod
# void volatility_gen::test (){
    def test():
        """

        Args:

        Returns:

        """
        #     string line ("#dynamic dbpedia 0.0 NORMAL frequency.txt volatility.txt 0.025");
        line = "#dynamic dbpedia 0.0 NORMAL frequency.txt volatility.txt 0.025"
#     //volatility_gen v_generator ("frequency.txt", 0.0, VOLATILITY_MODEL::CAUCHY_DIST, "volatility.txt");
#     string name;
#     float advance_pr;
#     volatility_gen * v_generator = volatility_gen::parse(line, name, advance_pr);
        v_generator = VolatilityGen.parse(line)  # , name, advance_pr)
#     cout << name << ", " << advance_pr << "\n";
#         print(name + ", " + advance_pr)

#     int model_size = 1000;
        model_size = 1000
#     v_generator->initialize(model_size);
        v_generator.initialize(model_size)

#     for (int i=0; i<10000; i++){
        for i in range(10000):
#         v_generator->advance ();
            v_generator.advance()
#     }

#     for (int j=0; j<100000; j++){
        for j in range(100000):
#         cout << v_generator->next_rand_index() << "\n";
            print(v_generator.next_rand_index())
#     }

#     for (int k=0; k<model_size; k++){
        for k in range(model_size):
#         cout << k << ":" << v_generator->get_probability(k) << "\n";
            print(f'{k}:{v_generator.get_probability(k)}')
#     }

#     delete v_generator;

#     /*
#     for (int i=0; i<10; i++){
#         cout << v_generator.get_volatility_value(i) << " ";
#         for (int j=0; j<1000; j++){
#             cout << v_generator.get_next_ln_increase(i) << " ";
#         }
#         cout << "\n";
#     }
#     */
# }
        pass  # end of volatility_gen::test

    pass  # end of VolatilityGen
