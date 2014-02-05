﻿/* Copyright (C) 2009-2014 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
                    string errMsg="Log-Mean Reverting calibration: input series cannot contain negative entries";
                   
                    return new EstimationResult(errMsg);
                }
            }

            //convert to Vector and transform to log-series
            Vector s = (Vector)Array.ConvertAll<IMarketData, double>(tmp, x => Math.Log( ((Scalar)x).Value));


            //Assumes elements are  order for decreasing date
            double dt = (tmp[0].TimeStamp - tmp[1].TimeStamp).TotalDays / 252.0;


            var calibrationResult=CalibrateOU(s, dt);

            
            // initial value and long term are expressed original value
            calibrationResult.Values[0] = Math.Exp(calibrationResult.Values[0]);
            calibrationResult.Values[2] = Math.Exp(calibrationResult.Values[2]);


            return calibrationResult;
        }

    }
}