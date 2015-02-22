/* Copyright (C) 2009-2014 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DVPLI;
using DVPLDOM;
using Mono.Addins;

namespace OrnsteinUhlenbeck
{
   
    /// <summary>
    /// Provides simulation for an Ornstein-Uhlenbeck process
    /// </summary>
    [Serializable]
    public class OrnsteinUhlenbeck : IExtensibleProcess, IFullSimulator, IParsable,IPopulable, IExportableContainer
    {
        static internal string processType = "Ornstein-Uhlenbeck";
        
        //Labels: do not change
        static string spotLabel = "Spot";
        static string muLabel = "Mu";
        static string sigmaLabel = "Sigma";
        static string lambdaLabel = "Lambda";



        public IModelParameter spot;
        public IModelParameter lambda;
        public IModelParameter sigma;
        public IModelParameter mu;

        /// <summary>
        /// Cached model parameters used for simulation
        /// </summary>
        [NonSerialized]
        double a, b, c;


        public OrnsteinUhlenbeck()
        {
            Init();
        }

        void Init()
        {
            this.spot = new ModelParameter(0,spotLabel);
            this.mu = new ModelParameter(0, muLabel);
            this.lambda = new ModelParameter(0, lambdaLabel);
            this.sigma = new ModelParameter(0, sigmaLabel);
        }


        public void Simulate(double[] Dates, IReadOnlyMatrixSlice Noise, IMatrixSlice OutDynamic)
        {
          OutDynamic[0,0]=spot.fV();
          for (int i = 1; i < Dates.Length; i++)
              OutDynamic[i, 0] = a * OutDynamic[i - 1, 0] + b + c * Noise[i-1, 0];
        }

        public bool ImplementsFullSimulation
        {
            get { return true; }
        }

        public bool ImplementsMarkovBasedSimulation
        {
            get { return false; }
        }

        public ProcessInfo ProcessInfo
        {
            get { return new ProcessInfo(processType); }
        }

        public bool Parse(IProject context)
        {
            bool errors = false;
            DVPLI.BoolHelper.AddBool(errors, spot.Parse(context));
            DVPLI.BoolHelper.AddBool(errors, lambda.Parse(context));
            DVPLI.BoolHelper.AddBool(errors, mu.Parse(context));
            DVPLI.BoolHelper.AddBool(errors, sigma.Parse(context));

            return errors;
        }

        public void Setup(double[] SimulationDates)
        {
            double dt = SimulationDates[1] - SimulationDates[0];//assumes constant dt

            a = Math.Exp(-lambda.fV() * dt);
            b = mu.fV() * (1 - Math.Exp(-lambda.fV() * dt));
            c = sigma.fV() * Math.Sqrt((1 - Math.Exp(-2 * lambda.fV() * dt)) / (2 * lambda.fV()));
        }

        public SimulationInfo SimulationInfo
        {
            get { var si=new SimulationInfo();
                    si.NoiseSize = 1;
                    si.StateSize = 1;
                    return si;
            }
        }

        public List<IExportable> ExportObjects(bool recursive)
        {
            List<IExportable> parameters = new List<IExportable>();
            parameters.Add(this.spot);
            parameters.Add(this.mu);
            parameters.Add(this.lambda);
            parameters.Add(this.sigma);
            return parameters;
        }


        public void Populate(string[] names, double[] values)
        {
            bool found;
            this.spot = new ModelParameter(PopulateHelper.GetValue(spotLabel, names, values, out found), spotLabel);
            this.mu = new ModelParameter(PopulateHelper.GetValue(muLabel,  names, values, out found), muLabel);
            this.lambda = new ModelParameter(PopulateHelper.GetValue(lambdaLabel, names, values, out found), lambdaLabel);
            this.sigma = new ModelParameter(PopulateHelper.GetValue(sigmaLabel,names, values, out found), sigmaLabel);
        }
    }
}
