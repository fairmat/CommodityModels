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
using Mono.Addins;
using DVPLI.MarketDataTypes;
using Fairmat.Calibration;
using Fairmat.MarketData;


namespace OrnsteinUhlenbeck
{


    /// <summary>
    /// Wrapper for the new model.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class OrnsteinUhlenbeckHistoricalCalibrationNew : OrnsteinUhlenbeckHistoricalCalibration
    {
        public override Type ProvidesTo
        {
            get { return typeof(OrnsteinUhlenbeck); }
        }
    }

    /// <summary>
    /// Calibrates the Ornstein Uhlenbeck from an historical price series
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class OrnsteinUhlenbeckHistoricalCalibration : IMenuItemDescription,IEstimatorEx
    {
        public IEstimationSettings DefaultSettings
        {
            get
            {
                return UserSettings.GetSettings(typeof(OrnsteinUhlenbeckCalibrationSettings)) as OrnsteinUhlenbeckCalibrationSettings;
            }
        }
        public virtual Type ProvidesTo
        {
            get { return typeof(DVPLDOM.StocasticProcessMeanReversion); }
        }

        public virtual string ToolTipText
        {
            get { return "Calibrate an Ornstein-Uhlenbeck from the stock/commodiy time series."; }
        }

        public virtual string Description
        {
            get { return "Ornstein-Uhlenbeck Historical Calibration"; }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="data"></param>
        /// <param name="settings"></param>
        /// <param name="controller"></param>
        /// <param name="properties"></param>
        /// <returns>An estimation results containing values for mu,sigma and lambda.</returns>
        public virtual EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            IMarketData[] tmp = data[0] as IMarketData[];
            tmp=tmp.OrderBy(x => x.TimeStamp).Reverse().ToArray();
           
            //convert to Vector

            Vector s = (Vector)Array.ConvertAll<IMarketData,double>(tmp, x => ((Scalar)x).Value);


            //Assume elements are  order for decreasing date
            double dt = (tmp[0].TimeStamp - tmp[1].TimeStamp).TotalDays / 252.0;


            return CalibrateOU(s, dt);
        }

        /// <summary>
        /// This code below replicate the calibration procedure reported below
        /// http://www.sitmo.com/article/calibrating-the-ornstein-uhlenbeck-model/
        /// </summary>
        /// <param name="s">The time series</param>
        /// <param name="delta">The average distance between two elements of the series.</param>
        /// <returns>The calibrated parameters</returns>
        protected static EstimationResult CalibrateOU(Vector s, double delta)
        {
            int n = s.Length - 1;
            double sx = 0;
            double sy = 0;
            double sxx = 0;
            double sxy = 0;
            double syy = 0;
            for (int i = 0; i < n; i++)
            {
                sx += s[i];
                sy += s[i + 1];
                sxx += s[i] * s[i];
                sxy += s[i] * s[i + 1];
                syy += s[i + 1] * s[i + 1];
            }
            //S_t=a S_{t-1}+ b+ sd epsilon
            double a = (n * sxy - sx * sy) / (n * sxx - Math.Pow(sx, 2));
            double b = (sy - a * sx) / n;
            double sd = Math.Sqrt((n * syy - Math.Pow(sy, 2) - a * (n * sxy - sx * sy)) / n / (n - 2));

            double spot = s[0];
            //Least Squares calibration
            double lambda = -Math.Log(a) / delta;
            double mu = b / (1 - a);
            double sigma = sd * Math.Sqrt(-2 * Math.Log(a) / delta / (1.0 - Math.Pow(a, 2)));


            //LL calibation
            mu = (sy * sxx - sx * sxy) / (n * (sxx - sxy) - (Math.Pow(sx, 2) - sx * sy));
            lambda = -Math.Log((sxy - mu * sx - mu * sy + n * Math.Pow(mu, 2)) / (sxx - 2 * mu * sx + n * Math.Pow(mu, 2))) / delta;
            a = Math.Exp(-lambda * delta);
            double sigmah2 = (syy - 2 * a * sxy + Math.Pow(a, 2) * sxx - 2 * mu * (1 - a) * (sy - a * sx) + n * Math.Pow(mu, 2) * Math.Pow(1 - a, 2)) / n;
            sigma = Math.Sqrt(sigmah2 * 2 * lambda / (1 - Math.Pow(a, 2)));

            /*
             * Parameters description
             * spot: the inital value
             * lambda: the mean reverting rate
             * mu:   the long term mean
             * sigma: the standard deviation
             */
            return new EstimationResult(new string[] { "spot", "lambda", "mu", "sigma" }, new double[] { spot, lambda, mu, sigma });
        }

        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[]{new EstimateRequirement(typeof(DVPLI.MarketDataTypes.Scalar[]))};
        }



     
    }
}
