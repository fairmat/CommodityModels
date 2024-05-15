﻿/* Copyright (C) 2009-2013 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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

using System.Reflection;
using Mono.Addins;

[assembly: AssemblyTrademark("")]
[assembly: AssemblyCulture("")]

// The following lines tell that the assembly is an addin
[assembly: Addin("Ornstein Uhlenbeck", "0.9.3", Category = "Calibration")]
[assembly: AddinDependency("Fairmat", "1.0")]
[assembly: AddinAuthor("Fairmat SRL")]
[assembly: AddinDescription("Ornstein Uhlenbeck / Vasicek / One Factor Mean Reverting  and Log-Mean Reverting models and historical calibration")]
