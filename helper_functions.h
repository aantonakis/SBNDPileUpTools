#ifndef PMT_HELPER_FUNCTIONS_H
#define PMT_HELPER_FUNCTIONS_H

// Standard Library Includes
#include <iostream>
#include <fstream>

// ROOT Includes
#include "TObject.h"
#include "TVector3.h"

namespace analysis {


// Function to find the indices of a vector associated with a particular value
std::vector<size_t> findIndices(const std::vector<int>& vec, int value) {
    std::vector<size_t> indices;
    auto it = vec.begin();

    while ((it = std::find(it, vec.end(), value)) != vec.end()) {
        indices.push_back(std::distance(vec.begin(), it));
        ++it; // Move iterator forward
    }

    return indices;
}


// function to calculate a vector of variances based on a window size
std::vector<double> rollingRMS(std::vector<unsigned int> waveform, long unsigned int w_size) {

    std::vector<double> rms_vals;
    long unsigned int i = 0;
    bool flag = true;
    while (flag) {
      if (i + w_size > waveform.size() - 1) {
        // reached the end of the waveform
        std::vector<unsigned int> w(i, waveform.size());
        double r = TMath::RMS(w.begin(), w.end());
        rms_vals.push_back(r);
        break;
      }
      std::vector<unsigned int> w(i, i + w_size);
      double r = TMath::RMS(w.begin(), w.end());
      rms_vals.push_back(r);
      i += w_size;
   }
   return rms_vals;
}


std::vector<double> discreteDerivative(std::vector<double> rms) {

    std::vector<double> d;
    for (long unsigned int i = 1; i < rms.size(); ++i) {
      d.push_back(rms.at(i) - rms.at(i-1));
    }
    return d;
}

// Potential Pulse Finding Algorithm
std::pair<std::vector<int>, std::vector<int>> find_pulses(std::vector<unsigned int> waveform, long unsigned int w_size) {

    std::vector<double> r = rollingRMS(waveform, w_size);
    std::vector<double> d = discreteDerivative(r);

    std::vector<int> starts;
    std::vector<int> ends;
    long unsigned int sample = 0;
    bool flag = false;
    while(sample < d.size()) {

      if (flag) {
        // look for the end of a pulse
        for (long unsigned int num = sample; num < d.size(); ++num) {
          int e = d.at(num);
          if ((-10 < e && e < 10) || (num + 1 >= d.size()) ) {

            ends.push_back(num*w_size);
            flag = false;
            sample = num + 1;
            break;
          }

        } // for loop to look for pulse end

      }
      // Found a pulse start?
      else if (d.at(sample) > 100 && !flag) {
        if (sample == d.size()-1) {
          sample += 1;
          break;
        }
        // Found the start of a pulse 
        int start = (sample + 1)*w_size;
        starts.push_back(start);
        flag = true;
        sample += 1;

      }
      else {
        //increment the sample
        sample += 1;

      }

    } // end of while loop

    std::pair<std::vector<int>, std::vector<int>> pulse(starts, ends);
    return pulse;

}


} // end of namespace

#endif
