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
using DVPLDOM;

namespace OrnsteinUhlenbeck
{
    //[Extension("/Fairmat/ProcessTypeChoice")] //Currently Disabled
    public class OrnsteinUhlenbeckChoice : IEditableChoice
    {
        #region IEditableOption Members

        /// <summary>
        /// Gets the name of the model which will be shown to the user.
        /// </summary>
        public string Description
        {
            get
            {
                return "Commodity/" + OrnsteinUhlenbeck.processType;
            }
        }

        /// <summary>
        /// Creates an IEditable instance from a StochasticProcessExtendible,
        /// which will handle the CIR plugin (<see cref="CIR"/> part).
        /// </summary>
        /// <returns>A reference to a new IEditable instance.</returns>
        public IEditable CreateInstance()
        {
            return new StochasticProcessExtendible(null, new OrnsteinUhlenbeck());
        }

        #endregion
    }


}
