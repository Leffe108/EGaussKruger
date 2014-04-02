EGaussKruger
============

This is a port of gausskruger.js to PHP packaged as a Yii extension. (but can easily be used for other PHP projects too)

gausskruger.js is a library to convert to/from common Swedish coordinate systems and the international WGS84 standard used by most GPS devices today. It can convert ``RT90 <-> WGS84`` and ``SWEREF99 <-> WGS84``.

The port is more or less just a translation from JavaScript to PHP with the code encapsulated inside a class. (the original library uses global variables, which is not the case here)

#Original JS version
Many thanks goes to the original JavaScript version that can be found here: 
http://mellifica.se/geodesi/gausskruger.js

#Usage
To convert from WGS84 to RT90 (2.5 gon west)
```php
$gk = new EGaussKruger();
$gk->swedish_params("rt90_2.5_gon_v");
list($rt90_x, $rt90_y) = $gk->geodetic_to_grid($latitude, $longitude);
```

To convert from RT90 (2.5 gon west) to WGS84
```php
$gk = new EGaussKruger();
$gk->swedish_params("rt90_2.5_gon_v");
list($latitude, $longitude) = $gk->grid_to_geodetic($rt90_x, $rt90_5);
```

If you open the EGaussKruger.php file and read the code for swedish_params you will see the exact name of each variant of RT90 and SWEREF99 that is supported.

#Use as Yii Extension
Just copy it to protected/extensions/ and make sure you autoload the protected/extensions/ directory.

#Use elsewhere
You can rename the class if you want to use it elsewhere where the E prefix is not wanted for extensions. It doesn't actually depend on anything in the Yii framework.
