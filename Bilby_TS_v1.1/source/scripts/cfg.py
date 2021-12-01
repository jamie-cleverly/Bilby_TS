#    cfg.py  Version management and specification
#    Copyright (C) 2016  Dr Jamie Cleverly
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
    Version History:
    <<1.1 6 June 2018, Update and new author name, Bilby TS Time Series Analyst>>
    <<1.0 6 June 2018, Update and new version name, Bilby TS Time Series Analyst>>
    
    Version History, OzFluxQC Simulator:
    <<2.9.5 20 December 2015, implementation of Chapin semantics for carbon fluxes (ammended)>>
    <<2.9.0 15 May 2015, Reconciled against OzFluxQC v2.9.0; variances, standard deviations and rotations updated>>
    <<2.8.6 2 Jan 2015, Reconciled against OzFluxQC v2.8.6>>
    <<v2.8.5 2 Oct 2014, Update following Alice Springs workshop>>
    <<v2.8.4a 18 Sep 2014, Reconciled against OzFluxQCv2.8.4, L1-L4 including new levels L4-L6>>
    <<v2.8.1a 12 Sep 2014, Reconciled against OzFluxQCv2.8.1, L1-L3 including new 3-D ncFile specifications>>
    <<v2.7 28 Oct 2013, Error in calculation of rhod corrected (rhod previously computed as total air density)>>
    <<v2.6 26 Sep 2013, Units error in WPLcov corrected>>
    <<v2.5 18 May 2013, Footprint model complete (rotated to geographic coordinates (m)-1d & 2d footprints, weighted footrpints, footprint climatology)>>
    <<v2.4 19 Apr 2013, WPL error introduced in v2.3 corrected, weighted footprint output added>>
    <<v2.3 24 Feb 2013, ET calculation added, controlfiles revised to simplify function calls>>
    <<v2.2 11 Feb 2013, Climatology implemented at L3 and L4; Penman-Monteith updated to output Gc, correct Gst computation to use air-to-soil q gradient rather than saturation deficit>>
    <<v2.1 3 Feb 2013, MOST footprint model (Kljun et al. 2004) implemented at L3>>
    <<v2.0 8 Jun 2012, version arrising from conclusion of OzFlux UTS data workshop>>
    <<v1.4 30 Sep 2011, final version arrising from OzFlux Black Mountain data workshop>>
    <<v1.0: 21 Jul 2011, code diversion reconciliation>>
"""

version_name = "Bilby_TS"
version_number = "v1.1"

# Bilby TS
# 1.1      - 1 December 2021, Update and new author name
# 1.0      - 6 June 2018, Update and new version name


# OzFluxQC Simulator
# 2.9.5    - 20 December 2015, implementation of Chapin semantics for carbon fluxes (ammended)
# 2.9.0    - 15 May 2015, Reconciled against OzFluxQC v2.9.0; variances, standard deviations and rotations updated
# 2.8.6    - 2 Jan 2015, Reconciled against OzFluxQC v2.8.6
# 2.8.5    - 2 Oct 2014, Update following Alice Springs workshop
# 2.8.4a   - 18 Sep 2014, Reconciled against OzFluxQCv2.8.4, L1-L4 including new levels L4-L6
# V2.8.1a  - 12 Sep 2013, Reconciled against OzFluxQCv2.8.1, L1-L3 including new 3-D ncFile specifications
# V2.7     - 28 Oct 2013, Error in calculation of rhod corrected
#            (rhod previously computed as total air density)
#          - cf conventions updated to match OzFluxQCv2.7.0
# V2.6     - 26 Sep 2013, Units error in WPLcov corrected
# V2.5     - 18 May 2013, Footprint model complete
#            (rotated to geographic coordinates (m)-1d & 2d footprints,
#            weighted footrpints, footprint climatology)
# V2.4     - 19 Apr 2013, WPL error introduced in v2.3 corrected
#          - weighted footprint output added
# V2.3     - 24 Feb 2013, ET calculation added
#          - controlfiles revised to simplify function calls
# V2.2     - 11 Feb 2013, Climatology implemented at L3 and L4
#          - Penman-Monteith updated to output Gc
#            correct Gst computation to use air-to-soil q gradient rather than saturation deficit
# V2.1     - 3 Feb 2013, MOST footprint model (Kljun et al. 2004) implemented at L3
# V2.0     - 8 Jun 2012, version arrising from conclusion of OzFlux UTS data workshop
# V1.4     - 30 Sep 2011, final version arrising from OzFlux Black Mountain data workshop
# V1.0     - 21 Jul 2011, code diversion reconciliation
