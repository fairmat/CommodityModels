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
using Mono.Addins;
using DVPLI;
using DVPLI.MarketDataTypes;
namespace OrnsteinUhlenbeck
{
    [Extension("/Fairmat/Estimator")]
    public class LogMeanRevertingHistoricalCalibration : OrnsteinUhlenbeckHistoricalCalibration
    {
        public override Type ProvidesTo
        {
            get { return typeof(DVPLDOM.StocasticProcessLogMeanReversion); }
        }

        public override EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            IMarketData[] tmp = data[0] as IMarketData[];
            tmp = tmp.OrderBy(x => x.TimeStamp).Reverse().ToArray();

            foreach (Scalar entry in tmp)
            {
                //check for positiveness
                if (entry.Value < 0)
                {
                    string errMsg = "Log-Mean Reverting calibration: input series cannot contain negative entries";

                    return new EstimationResult(errMsg);
                }
            }

            //convert to Vector and transform to log-series
            Vector s = (Vector)Array.ConvertAll<IMarketData, double>(tmp, x => Math.Log(((Scalar)x).Value));


            //Assumes elements are  order for decreasing date
            double dt = (tmp[0].TimeStamp - tmp[1].TimeStamp).TotalDays / 252.0;


            var calibrationResult = CalibrateOU(s, dt);


            // initial value and long term are expressed original value
            calibrationResult.Values[0] = Math.Exp(calibrationResult.Values[0]);
            calibrationResult.Values[2] = Math.Exp(calibrationResult.Values[2]);


            return calibrationResult;
        }

    }

    /// <summary>
    /// This calibration procedure calibrates volatility and mean reversion to historical values
    /// but set long-term value and initial value to a fraction (longTermDeprecation) of the current value.
    /// </summary>
    [Obsolete]
    [Extension("/Fairmat/Estimator")]
    public class LongTermLogMeanRevertingHistoricalCalibration : LogMeanRevertingHistoricalCalibration
    {
        static double longTermDeprecation = 0.8;

        public override EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            // Use historical calibration
            var calibrationResult = base.Estimate(data, settings, controller, properties);

            // Set long-term and initial value to the same value
            calibrationResult.Values[0] = longTermDeprecation * calibrationResult.Values[0];
            calibrationResult.Values[2] = calibrationResult.Values[0];

            return calibrationResult;
        }

        public override string Description
        {
            get { return "Long-Term Ornstein-Uhlenbeck Historical Calibration"; }
        }
    }

    /// <summary>
    /// Specializes the calibration process by forcing a market data request using the DividendYield field.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class LogMeanRevertingHistoricalCalibrationForDividendYields : LogMeanRevertingHistoricalCalibration
    {
        public override string Description
        {
            get { return "Dividend Yield - Ornstein-Uhlenbeck Historical Calibration"; }
        }

        public override EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(DVPLI.MarketDataTypes.Scalar[])) { Field="DividendYield"} };
        }


    }
}



