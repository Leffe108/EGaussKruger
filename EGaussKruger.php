<?php 
/* -- Original copyright from http://mellifica.se/geodesi/gausskruger.js ---- */

// 
// Author: Arnold Andreasson, info@mellifica.se
// Copyright (c) 2007-2013 Arnold Andreasson 
// License: MIT License as follows:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// =============================================================================
// Javascript-implementation of "Gauss Conformal Projection 
// (Transverse Mercator), Krügers Formulas".
// - Parameters for SWEREF99 lat-long to/from RT90 and SWEREF99 
//   coordinates (RT90 and SWEREF99 are used in Swedish maps).
// Source: http://www.lantmateriet.se/geodesi/

/* -- This is EGaussKruger - Yii extension port -- */
//
// Port author: Leif Linse (Leffe108 on github)
// License: MIT License (see above)
//
// The port is more or less just a translation from JavaScript 
// into PHP with the code being encapsulated into a class.
// Many thanks to the original author who has made a much
// larger job in making the initial JavaScript version, than
// this quite simple port.

class EGaussKruger
{
	protected $axis = null; // Semi-major axis of the ellipsoid.
	protected $flattening = null; // Flattening of the ellipsoid.
	protected $central_meridian = null; // Central meridian for the projection.
	protected $lat_of_origin = null; // Latitude of origin.
	protected $scale = null; // Scale on central meridian.
	protected $false_northing = null; // Offset for origo.
	protected $false_easting = null; // Offset for origo.

	// Parameters for RT90 and SWEREF99TM.
	// Note: Parameters for RT90 are choosen to eliminate the 
	// differences between Bessel and GRS80-ellipsoides.
	// Bessel-variants should only be used if lat/long are given as
	// RT90-lat/long based on the Bessel ellipsoide (from old maps).
	// Parameter: projection (string). Must match if-statement.
	public function swedish_params($projection) {
		// RT90 parameters, GRS 80 ellipsoid.
		if ($projection == "rt90_7.5_gon_v") {
			$this->grs80_params();
			$this->central_meridian = 11.0 + 18.375/60.0;
			$this->scale = 1.000006000000;
			$this->false_northing = -667.282;
			$this->false_easting = 1500025.141;
		}
		else if ($projection == "rt90_5.0_gon_v") {
			$this->grs80_params();
			$this->central_meridian = 13.0 + 33.376/60.0;
			$this->scale = 1.000005800000;
			$this->false_northing = -667.130;
			$this->false_easting = 1500044.695;
		}
		else if ($projection == "rt90_2.5_gon_v") {
			$this->grs80_params();
			$this->central_meridian = 15.0 + 48.0/60.0 + 22.624306/3600.0;
			$this->scale = 1.00000561024;
			$this->false_northing = -667.711;
			$this->false_easting = 1500064.274;
		}
		else if ($projection == "rt90_0.0_gon_v") {
			$this->grs80_params();
			$this->central_meridian = 18.0 + 3.378/60.0;
			$this->scale = 1.000005400000;
			$this->false_northing = -668.844;
			$this->false_easting = 1500083.521;
		}
		else if ($projection == "rt90_2.5_gon_o") {
			$this->grs80_params();
			$this->central_meridian = 20.0 + 18.379/60.0;
			$this->scale = 1.000005200000;
			$this->false_northing = -670.706;
			$this->false_easting = 1500102.765;
		}
		else if ($projection == "rt90_5.0_gon_o") {
			$this->grs80_params();
			$this->central_meridian = 22.0 + 33.380/60.0;
			$this->scale = 1.000004900000;
			$this->false_northing = -672.557;
			$this->false_easting = 1500121.846;
		}
		
		// RT90 parameters, Bessel 1841 ellipsoid.
		else if ($projection == "bessel_rt90_7.5_gon_v") {
			$this->bessel_params();
			$this->central_meridian = 11.0 + 18.0/60.0 + 29.8/3600.0;
		}
		else if ($projection == "bessel_rt90_5.0_gon_v") {
			$this->bessel_params();
			$this->central_meridian = 13.0 + 33.0/60.0 + 29.8/3600.0;
		}
		else if ($projection == "bessel_rt90_2.5_gon_v") {
			$this->bessel_params();
			$this->central_meridian = 15.0 + 48.0/60.0 + 29.8/3600.0;
		}
		else if ($projection == "bessel_rt90_0.0_gon_v") {
			$this->bessel_params();
			$this->central_meridian = 18.0 + 3.0/60.0 + 29.8/3600.0;
		}
		else if ($projection == "bessel_rt90_2.5_gon_o") {
			$this->bessel_params();
			$this->central_meridian = 20.0 + 18.0/60.0 + 29.8/3600.0;
		}
		else if ($projection == "bessel_rt90_5.0_gon_o") {
			$this->bessel_params();
			$this->central_meridian = 22.0 + 33.0/60.0 + 29.8/3600.0;
		}

		// SWEREF99TM and SWEREF99ddmm  parameters.
		else if ($projection == "sweref_99_tm") {
			$this->sweref99_params();
			$this->central_meridian = 15.00;
			$this->lat_of_origin = 0.0;
			$this->scale = 0.9996;
			$this->false_northing = 0.0;
			$this->false_easting = 500000.0;
		}
		else if ($projection == "sweref_99_1200") {
			$this->sweref99_params();
			$this->central_meridian = 12.00;
		}
		else if ($projection == "sweref_99_1330") {
			$this->sweref99_params();
			$this->central_meridian = 13.50;
		}
		else if ($projection == "sweref_99_1500") {
			$this->sweref99_params();
			$this->central_meridian = 15.00;
		}
		else if ($projection == "sweref_99_1630") {
			$this->sweref99_params();
			$this->central_meridian = 16.50;
		}
		else if ($projection == "sweref_99_1800") {
			$this->sweref99_params();
			$this->central_meridian = 18.00;
		}
		else if ($projection == "sweref_99_1415") {
			$this->sweref99_params();
			$this->central_meridian = 14.25;
		}
		else if ($projection == "sweref_99_1545") {
			$this->sweref99_params();
			$this->central_meridian = 15.75;
		}
		else if ($projection == "sweref_99_1715") {
			$this->sweref99_params();
			$this->central_meridian = 17.25;
		}
		else if ($projection == "sweref_99_1845") {
			$this->sweref99_params();
			$this->central_meridian = 18.75;
		}
		else if ($projection == "sweref_99_2015") {
			$this->sweref99_params();
			$this->central_meridian = 20.25;
		}
		else if ($projection == "sweref_99_2145") {
			$this->sweref99_params();
			$this->central_meridian = 21.75;
		}
		else if ($projection == "sweref_99_2315") {
			$this->sweref99_params();
			$this->central_meridian = 23.25;
		}

		// Test-case:
		//	Lat: 66 0'0", lon: 24 0'0".
		//	X:1135809.413803 Y:555304.016555.
		else if ($projection == "test_case") {
			$this->axis = 6378137.0;
			$this->flattening = 1.0 / 298.257222101;
			$this->central_meridian = 13.0 + 35.0/60.0 + 7.692000/3600.0;
			$this->lat_of_origin = 0.0;
			$this->scale = 1.000002540000;
			$this->false_northing = -6226307.8640;
			$this->false_easting = 84182.8790;

		// Not a valid projection.		
		} else {
			$this->central_meridian = null;
		}
	}
	// Sets of default parameters.
	protected function grs80_params() {
		$this->axis = 6378137.0; // GRS 80.
		$this->flattening = 1.0 / 298.257222101; // GRS 80.
		$this->central_meridian = null;
		$this->lat_of_origin = 0.0;
	}
	protected function bessel_params() {
		$this->axis = 6377397.155; // Bessel 1841.
		$this->flattening = 1.0 / 299.1528128; // Bessel 1841.
		$this->central_meridian = null;
		$this->lat_of_origin = 0.0;
		$this->scale = 1.0;
		$this->false_northing = 0.0;
		$this->false_easting = 1500000.0;
	}
	protected function sweref99_params() {
		$this->axis = 6378137.0; // GRS 80.
		$this->flattening = 1.0 / 298.257222101; // GRS 80.
		$this->central_meridian = null;
		$this->lat_of_origin = 0.0;
		$this->scale = 1.0;
		$this->false_northing = 0.0;
		$this->false_easting = 150000.0;
	}

	// Conversion from geodetic coordinates to grid coordinates.
	public function geodetic_to_grid($latitude, $longitude) {
		$x_y = array(null, null);
		if ($this->central_meridian == null) {
			return $x_y;
		}
		// Prepare ellipsoid-based stuff.
		$e2 = $this->flattening * (2.0 - $this->flattening);
		$n = $this->flattening / (2.0 - $this->flattening);
		$a_roof = $this->axis / (1.0 + $n) * (1.0 + $n*$n/4.0 + $n*$n*$n*$n/64.0);
		$A = $e2;
		$B = (5.0*$e2*$e2 - $e2*$e2*$e2) / 6.0;
		$C = (104.0*$e2*$e2*$e2 - 45.0*$e2*$e2*$e2*$e2) / 120.0;
		$D = (1237.0*$e2*$e2*$e2*$e2) / 1260.0;
		$beta1 = $n/2.0 - 2.0*$n*$n/3.0 + 5.0*$n*$n*$n/16.0 + 41.0*$n*$n*$n*$n/180.0;
		$beta2 = 13.0*$n*$n/48.0 - 3.0*$n*$n*$n/5.0 + 557.0*$n*$n*$n*$n/1440.0;
		$beta3 = 61.0*$n*$n*$n/240.0 - 103.0*$n*$n*$n*$n/140.0;
		$beta4 = 49561.0*$n*$n*$n*$n/161280.0;
		
		// Convert.
		$deg_to_rad = pi() / 180.0;
		$phi = $latitude * $deg_to_rad;
		$lambda = $longitude * $deg_to_rad;
		$lambda_zero = $this->central_meridian * $deg_to_rad;
		
		$phi_star = $phi - sin($phi) * cos($phi) * ($A + 
						$B*pow(sin($phi), 2) + 
						$C*pow(sin($phi), 4) + 
						$D*pow(sin($phi), 6));
		$delta_lambda = $lambda - $lambda_zero;
		$xi_prim = atan(tan($phi_star) / cos($delta_lambda));
		$eta_prim = atanh(cos($phi_star) * sin($delta_lambda));
		$x = $this->scale * $a_roof * ($xi_prim +
						$beta1 * sin(2.0*$xi_prim) * cosh(2.0*$eta_prim) +
						$beta2 * sin(4.0*$xi_prim) * cosh(4.0*$eta_prim) +
						$beta3 * sin(6.0*$xi_prim) * cosh(6.0*$eta_prim) +
						$beta4 * sin(8.0*$xi_prim) * cosh(8.0*$eta_prim)) + 
						$this->false_northing;
		$y = $this->scale * $a_roof * ($eta_prim +
						$beta1 * cos(2.0*$xi_prim) * sinh(2.0*$eta_prim) +
						$beta2 * cos(4.0*$xi_prim) * sinh(4.0*$eta_prim) +
						$beta3 * cos(6.0*$xi_prim) * sinh(6.0*$eta_prim) +
						$beta4 * cos(8.0*$xi_prim) * sinh(8.0*$eta_prim)) + 
						$this->false_easting;
		$x_y[0] = round($x * 1000.0) / 1000.0;
		$x_y[1] = round($y * 1000.0) / 1000.0;
	//	x_y[0] = x;
	//	x_y[1] = y;
		return $x_y;
	}

	// Conversion from grid coordinates to geodetic coordinates.
	public function grid_to_geodetic($x, $y) {
		$lat_lon = array(null, null);
		if ($this->central_meridian == null) {
			return $lat_lon;
		}
		// Prepare ellipsoid-based stuff.
		$e2 = $this->flattening * (2.0 - $this->flattening);
		$n = $this->flattening / (2.0 - $this->flattening);
		$a_roof = $this->axis / (1.0 + $n) * (1.0 + $n*$n/4.0 + $n*$n*$n*$n/64.0);
		$delta1 = $n/2.0 - 2.0*$n*$n/3.0 + 37.0*$n*$n*$n/96.0 - $n*$n*$n*$n/360.0;
		$delta2 = $n*$n/48.0 + $n*$n*$n/15.0 - 437.0*$n*$n*$n*$n/1440.0;
		$delta3 = 17.0*$n*$n*$n/480.0 - 37*$n*$n*$n*$n/840.0;
		$delta4 = 4397.0*$n*$n*$n*$n/161280.0;
		
		$Astar = $e2 + $e2*$e2 + $e2*$e2*$e2 + $e2*$e2*$e2*$e2;
		$Bstar = -(7.0*$e2*$e2 + 17.0*$e2*$e2*$e2 + 30.0*$e2*$e2*$e2*$e2) / 6.0;
		$Cstar = (224.0*$e2*$e2*$e2 + 889.0*$e2*$e2*$e2*$e2) / 120.0;
		$Dstar = -(4279.0*$e2*$e2*$e2*$e2) / 1260.0;

		// Convert.
		$deg_to_rad = pi() / 180;
		$lambda_zero = $this->central_meridian * $deg_to_rad;
		$xi = ($x - $this->false_northing) / ($this->scale * $a_roof);		
		$eta = ($y - $this->false_easting) / ($this->scale * $a_roof);
		$xi_prim = $xi - 
						$delta1*sin(2.0*$xi) * cosh(2.0*$eta) - 
						$delta2*sin(4.0*$xi) * cosh(4.0*$eta) - 
						$delta3*sin(6.0*$xi) * cosh(6.0*$eta) - 
						$delta4*sin(8.0*$xi) * cosh(8.0*$eta);
		$eta_prim = $eta - 
						$delta1*cos(2.0*$xi) * sinh(2.0*$eta) - 
						$delta2*cos(4.0*$xi) * sinh(4.0*$eta) - 
						$delta3*cos(6.0*$xi) * sinh(6.0*$eta) - 
						$delta4*cos(8.0*$xi) * sinh(8.0*$eta);
		$phi_star = asin(sin($xi_prim) / cosh($eta_prim));
		$delta_lambda = atan(sinh($eta_prim) / cos($xi_prim));
		$lon_radian = $lambda_zero + $delta_lambda;
		$lat_radian = $phi_star + sin($phi_star) * cos($phi_star) * 
						($Astar + 
						 $Bstar*pow(sin($phi_star), 2) + 
						 $Cstar*pow(sin($phi_star), 4) + 
						 $Dstar*pow(sin($phi_star), 6));  	
		$lat_lon[0] = $lat_radian * 180.0 / pi();
		$lat_lon[1] = $lon_radian * 180.0 / pi();
		return $lat_lon;
	}
}
