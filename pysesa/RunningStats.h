/* PySESA (Python program for Spatially Explicit Spectral Anlysis) 
 has been developed at the Grand Canyon Monitorinf & Research Center,
 U.S. Geological Survey

 Author: Daniel Buscombe
 Project homepage: <https://github.com/dbuscombe-usgs/pysesa>

This software is in the public domain because it contains materials that originally came from 
the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

   http://www.johndcook.com/blog/skewness_kurtosis/ */

#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H
 
class RunningStats
{
public:
    RunningStats();
    void Clear();
    void Push(double x);
    long long NumDataValues() const;
    double Mean() const;
    double Variance() const;
    double StandardDeviation() const;
    double Skewness() const;
    double Kurtosis() const;
 
    friend RunningStats operator+(const RunningStats a, const RunningStats b);
    RunningStats& operator+=(const RunningStats &rhs);
 
private:
    long long n;
    double M1, M2, M3, M4;
};
 
#endif
