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
using Fairmat.MarketData;
namespace OrnsteinUhlenbeck
{
    [SettingsContainer("Ornstein Uhlenbeck Calibration Preferences", SettingType = SettingType.Calibration)]
    [Mono.Addins.Extension("/Fairmat/UserSettings")]
    [Serializable]
    public class OrnsteinUhlenbeckCalibrationSettings : TimeSeriesEstimationSettings
    {
        public OrnsteinUhlenbeckCalibrationSettings()
        {
            this.Freq = SampleFrequency.Daily;
            this.MinimumNumberOfSamples = 1;
            this.MaxRequestedSamples = 252;
        }
    }


}
