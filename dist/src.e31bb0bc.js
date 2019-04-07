// modules are defined as an array
// [ module function, map of requires ]
//
// map of requires is short require name -> numeric require
//
// anything defined in a previous bundle is accessed via the
// orig method which is the require for previous bundles
parcelRequire = (function (modules, cache, entry, globalName) {
  // Save the require from previous bundle to this closure if any
  var previousRequire = typeof parcelRequire === 'function' && parcelRequire;
  var nodeRequire = typeof require === 'function' && require;

  function newRequire(name, jumped) {
    if (!cache[name]) {
      if (!modules[name]) {
        // if we cannot find the module within our internal map or
        // cache jump to the current global require ie. the last bundle
        // that was added to the page.
        var currentRequire = typeof parcelRequire === 'function' && parcelRequire;
        if (!jumped && currentRequire) {
          return currentRequire(name, true);
        }

        // If there are other bundles on this page the require from the
        // previous one is saved to 'previousRequire'. Repeat this as
        // many times as there are bundles until the module is found or
        // we exhaust the require chain.
        if (previousRequire) {
          return previousRequire(name, true);
        }

        // Try the node require function if it exists.
        if (nodeRequire && typeof name === 'string') {
          return nodeRequire(name);
        }

        var err = new Error('Cannot find module \'' + name + '\'');
        err.code = 'MODULE_NOT_FOUND';
        throw err;
      }

      localRequire.resolve = resolve;
      localRequire.cache = {};

      var module = cache[name] = new newRequire.Module(name);

      modules[name][0].call(module.exports, localRequire, module, module.exports, this);
    }

    return cache[name].exports;

    function localRequire(x){
      return newRequire(localRequire.resolve(x));
    }

    function resolve(x){
      return modules[name][1][x] || x;
    }
  }

  function Module(moduleName) {
    this.id = moduleName;
    this.bundle = newRequire;
    this.exports = {};
  }

  newRequire.isParcelRequire = true;
  newRequire.Module = Module;
  newRequire.modules = modules;
  newRequire.cache = cache;
  newRequire.parent = previousRequire;
  newRequire.register = function (id, exports) {
    modules[id] = [function (require, module) {
      module.exports = exports;
    }, {}];
  };

  var error;
  for (var i = 0; i < entry.length; i++) {
    try {
      newRequire(entry[i]);
    } catch (e) {
      // Save first error but execute all entries
      if (!error) {
        error = e;
      }
    }
  }

  if (entry.length) {
    // Expose entry point to Node, AMD or browser globals
    // Based on https://github.com/ForbesLindesay/umd/blob/master/template.js
    var mainExports = newRequire(entry[entry.length - 1]);

    // CommonJS
    if (typeof exports === "object" && typeof module !== "undefined") {
      module.exports = mainExports;

    // RequireJS
    } else if (typeof define === "function" && define.amd) {
     define(function () {
       return mainExports;
     });

    // <script>
    } else if (globalName) {
      this[globalName] = mainExports;
    }
  }

  // Override the current require with this new one
  parcelRequire = newRequire;

  if (error) {
    // throw error from earlier, _after updating parcelRequire_
    throw error;
  }

  return newRequire;
})({"../node_modules/d3-geo/src/adder.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

// Adds floating point numbers with twice the normal precision.
// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
// 305–363 (1997).
// Code adapted from GeographicLib by Charles F. F. Karney,
// http://geographiclib.sourceforge.net/
function _default() {
  return new Adder();
}

function Adder() {
  this.reset();
}

Adder.prototype = {
  constructor: Adder,
  reset: function () {
    this.s = // rounded value
    this.t = 0; // exact error
  },
  add: function (y) {
    add(temp, y, this.t);
    add(this, temp.s, this.s);
    if (this.s) this.t += temp.t;else this.s = temp.t;
  },
  valueOf: function () {
    return this.s;
  }
};
var temp = new Adder();

function add(adder, a, b) {
  var x = adder.s = a + b,
      bv = x - a,
      av = x - bv;
  adder.t = a - av + (b - bv);
}
},{}],"../node_modules/d3-geo/src/math.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.acos = acos;
exports.asin = asin;
exports.haversin = haversin;
exports.tan = exports.sqrt = exports.sign = exports.sin = exports.pow = exports.log = exports.floor = exports.exp = exports.ceil = exports.cos = exports.atan2 = exports.atan = exports.abs = exports.radians = exports.degrees = exports.tau = exports.quarterPi = exports.halfPi = exports.pi = exports.epsilon2 = exports.epsilon = void 0;
var epsilon = 1e-6;
exports.epsilon = epsilon;
var epsilon2 = 1e-12;
exports.epsilon2 = epsilon2;
var pi = Math.PI;
exports.pi = pi;
var halfPi = pi / 2;
exports.halfPi = halfPi;
var quarterPi = pi / 4;
exports.quarterPi = quarterPi;
var tau = pi * 2;
exports.tau = tau;
var degrees = 180 / pi;
exports.degrees = degrees;
var radians = pi / 180;
exports.radians = radians;
var abs = Math.abs;
exports.abs = abs;
var atan = Math.atan;
exports.atan = atan;
var atan2 = Math.atan2;
exports.atan2 = atan2;
var cos = Math.cos;
exports.cos = cos;
var ceil = Math.ceil;
exports.ceil = ceil;
var exp = Math.exp;
exports.exp = exp;
var floor = Math.floor;
exports.floor = floor;
var log = Math.log;
exports.log = log;
var pow = Math.pow;
exports.pow = pow;
var sin = Math.sin;
exports.sin = sin;

var sign = Math.sign || function (x) {
  return x > 0 ? 1 : x < 0 ? -1 : 0;
};

exports.sign = sign;
var sqrt = Math.sqrt;
exports.sqrt = sqrt;
var tan = Math.tan;
exports.tan = tan;

function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}

function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}

function haversin(x) {
  return (x = sin(x / 2)) * x;
}
},{}],"../node_modules/d3-geo/src/noop.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = noop;

function noop() {}
},{}],"../node_modules/d3-geo/src/stream.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function streamGeometry(geometry, stream) {
  if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
    streamGeometryType[geometry.type](geometry, stream);
  }
}

var streamObjectType = {
  Feature: function (object, stream) {
    streamGeometry(object.geometry, stream);
  },
  FeatureCollection: function (object, stream) {
    var features = object.features,
        i = -1,
        n = features.length;

    while (++i < n) streamGeometry(features[i].geometry, stream);
  }
};
var streamGeometryType = {
  Sphere: function (object, stream) {
    stream.sphere();
  },
  Point: function (object, stream) {
    object = object.coordinates;
    stream.point(object[0], object[1], object[2]);
  },
  MultiPoint: function (object, stream) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) object = coordinates[i], stream.point(object[0], object[1], object[2]);
  },
  LineString: function (object, stream) {
    streamLine(object.coordinates, stream, 0);
  },
  MultiLineString: function (object, stream) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) streamLine(coordinates[i], stream, 0);
  },
  Polygon: function (object, stream) {
    streamPolygon(object.coordinates, stream);
  },
  MultiPolygon: function (object, stream) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) streamPolygon(coordinates[i], stream);
  },
  GeometryCollection: function (object, stream) {
    var geometries = object.geometries,
        i = -1,
        n = geometries.length;

    while (++i < n) streamGeometry(geometries[i], stream);
  }
};

function streamLine(coordinates, stream, closed) {
  var i = -1,
      n = coordinates.length - closed,
      coordinate;
  stream.lineStart();

  while (++i < n) coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);

  stream.lineEnd();
}

function streamPolygon(coordinates, stream) {
  var i = -1,
      n = coordinates.length;
  stream.polygonStart();

  while (++i < n) streamLine(coordinates[i], stream, 1);

  stream.polygonEnd();
}

function _default(object, stream) {
  if (object && streamObjectType.hasOwnProperty(object.type)) {
    streamObjectType[object.type](object, stream);
  } else {
    streamGeometry(object, stream);
  }
}
},{}],"../node_modules/d3-geo/src/area.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.areaStream = exports.areaRingSum = void 0;

var _adder = _interopRequireDefault(require("./adder"));

var _math = require("./math");

var _noop = _interopRequireDefault(require("./noop"));

var _stream = _interopRequireDefault(require("./stream"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var areaRingSum = (0, _adder.default)();
exports.areaRingSum = areaRingSum;
var areaSum = (0, _adder.default)(),
    lambda00,
    phi00,
    lambda0,
    cosPhi0,
    sinPhi0;
var areaStream = {
  point: _noop.default,
  lineStart: _noop.default,
  lineEnd: _noop.default,
  polygonStart: function () {
    areaRingSum.reset();
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function () {
    var areaRing = +areaRingSum;
    areaSum.add(areaRing < 0 ? _math.tau + areaRing : areaRing);
    this.lineStart = this.lineEnd = this.point = _noop.default;
  },
  sphere: function () {
    areaSum.add(_math.tau);
  }
};
exports.areaStream = areaStream;

function areaRingStart() {
  areaStream.point = areaPointFirst;
}

function areaRingEnd() {
  areaPoint(lambda00, phi00);
}

function areaPointFirst(lambda, phi) {
  areaStream.point = areaPoint;
  lambda00 = lambda, phi00 = phi;
  lambda *= _math.radians, phi *= _math.radians;
  lambda0 = lambda, cosPhi0 = (0, _math.cos)(phi = phi / 2 + _math.quarterPi), sinPhi0 = (0, _math.sin)(phi);
}

function areaPoint(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  phi = phi / 2 + _math.quarterPi; // half the angular distance from south pole
  // Spherical excess E for a spherical triangle with vertices: south pole,
  // previous point, current point.  Uses a formula derived from Cagnoli’s
  // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).

  var dLambda = lambda - lambda0,
      sdLambda = dLambda >= 0 ? 1 : -1,
      adLambda = sdLambda * dLambda,
      cosPhi = (0, _math.cos)(phi),
      sinPhi = (0, _math.sin)(phi),
      k = sinPhi0 * sinPhi,
      u = cosPhi0 * cosPhi + k * (0, _math.cos)(adLambda),
      v = k * sdLambda * (0, _math.sin)(adLambda);
  areaRingSum.add((0, _math.atan2)(v, u)); // Advance the previous points.

  lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
}

function _default(object) {
  areaSum.reset();
  (0, _stream.default)(object, areaStream);
  return areaSum * 2;
}
},{"./adder":"../node_modules/d3-geo/src/adder.js","./math":"../node_modules/d3-geo/src/math.js","./noop":"../node_modules/d3-geo/src/noop.js","./stream":"../node_modules/d3-geo/src/stream.js"}],"../node_modules/d3-geo/src/cartesian.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.spherical = spherical;
exports.cartesian = cartesian;
exports.cartesianDot = cartesianDot;
exports.cartesianCross = cartesianCross;
exports.cartesianAddInPlace = cartesianAddInPlace;
exports.cartesianScale = cartesianScale;
exports.cartesianNormalizeInPlace = cartesianNormalizeInPlace;

var _math = require("./math");

function spherical(cartesian) {
  return [(0, _math.atan2)(cartesian[1], cartesian[0]), (0, _math.asin)(cartesian[2])];
}

function cartesian(spherical) {
  var lambda = spherical[0],
      phi = spherical[1],
      cosPhi = (0, _math.cos)(phi);
  return [cosPhi * (0, _math.cos)(lambda), cosPhi * (0, _math.sin)(lambda), (0, _math.sin)(phi)];
}

function cartesianDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cartesianCross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
} // TODO return a


function cartesianAddInPlace(a, b) {
  a[0] += b[0], a[1] += b[1], a[2] += b[2];
}

function cartesianScale(vector, k) {
  return [vector[0] * k, vector[1] * k, vector[2] * k];
} // TODO return d


function cartesianNormalizeInPlace(d) {
  var l = (0, _math.sqrt)(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  d[0] /= l, d[1] /= l, d[2] /= l;
}
},{"./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/bounds.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _adder = _interopRequireDefault(require("./adder"));

var _area = require("./area");

var _cartesian = require("./cartesian");

var _math = require("./math");

var _stream = _interopRequireDefault(require("./stream"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var lambda0,
    phi0,
    lambda1,
    phi1,
    // bounds
lambda2,
    // previous lambda-coordinate
lambda00,
    phi00,
    // first point
p0,
    // previous 3D point
deltaSum = (0, _adder.default)(),
    ranges,
    range;
var boundsStream = {
  point: boundsPoint,
  lineStart: boundsLineStart,
  lineEnd: boundsLineEnd,
  polygonStart: function () {
    boundsStream.point = boundsRingPoint;
    boundsStream.lineStart = boundsRingStart;
    boundsStream.lineEnd = boundsRingEnd;
    deltaSum.reset();

    _area.areaStream.polygonStart();
  },
  polygonEnd: function () {
    _area.areaStream.polygonEnd();

    boundsStream.point = boundsPoint;
    boundsStream.lineStart = boundsLineStart;
    boundsStream.lineEnd = boundsLineEnd;
    if (_area.areaRingSum < 0) lambda0 = -(lambda1 = 180), phi0 = -(phi1 = 90);else if (deltaSum > _math.epsilon) phi1 = 90;else if (deltaSum < -_math.epsilon) phi0 = -90;
    range[0] = lambda0, range[1] = lambda1;
  }
};

function boundsPoint(lambda, phi) {
  ranges.push(range = [lambda0 = lambda, lambda1 = lambda]);
  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
}

function linePoint(lambda, phi) {
  var p = (0, _cartesian.cartesian)([lambda * _math.radians, phi * _math.radians]);

  if (p0) {
    var normal = (0, _cartesian.cartesianCross)(p0, p),
        equatorial = [normal[1], -normal[0], 0],
        inflection = (0, _cartesian.cartesianCross)(equatorial, normal);
    (0, _cartesian.cartesianNormalizeInPlace)(inflection);
    inflection = (0, _cartesian.spherical)(inflection);
    var delta = lambda - lambda2,
        sign = delta > 0 ? 1 : -1,
        lambdai = inflection[0] * _math.degrees * sign,
        phii,
        antimeridian = (0, _math.abs)(delta) > 180;

    if (antimeridian ^ (sign * lambda2 < lambdai && lambdai < sign * lambda)) {
      phii = inflection[1] * _math.degrees;
      if (phii > phi1) phi1 = phii;
    } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign * lambda2 < lambdai && lambdai < sign * lambda)) {
      phii = -inflection[1] * _math.degrees;
      if (phii < phi0) phi0 = phii;
    } else {
      if (phi < phi0) phi0 = phi;
      if (phi > phi1) phi1 = phi;
    }

    if (antimeridian) {
      if (lambda < lambda2) {
        if (angle(lambda0, lambda) > angle(lambda0, lambda1)) lambda1 = lambda;
      } else {
        if (angle(lambda, lambda1) > angle(lambda0, lambda1)) lambda0 = lambda;
      }
    } else {
      if (lambda1 >= lambda0) {
        if (lambda < lambda0) lambda0 = lambda;
        if (lambda > lambda1) lambda1 = lambda;
      } else {
        if (lambda > lambda2) {
          if (angle(lambda0, lambda) > angle(lambda0, lambda1)) lambda1 = lambda;
        } else {
          if (angle(lambda, lambda1) > angle(lambda0, lambda1)) lambda0 = lambda;
        }
      }
    }
  } else {
    ranges.push(range = [lambda0 = lambda, lambda1 = lambda]);
  }

  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
  p0 = p, lambda2 = lambda;
}

function boundsLineStart() {
  boundsStream.point = linePoint;
}

function boundsLineEnd() {
  range[0] = lambda0, range[1] = lambda1;
  boundsStream.point = boundsPoint;
  p0 = null;
}

function boundsRingPoint(lambda, phi) {
  if (p0) {
    var delta = lambda - lambda2;
    deltaSum.add((0, _math.abs)(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
  } else {
    lambda00 = lambda, phi00 = phi;
  }

  _area.areaStream.point(lambda, phi);

  linePoint(lambda, phi);
}

function boundsRingStart() {
  _area.areaStream.lineStart();
}

function boundsRingEnd() {
  boundsRingPoint(lambda00, phi00);

  _area.areaStream.lineEnd();

  if ((0, _math.abs)(deltaSum) > _math.epsilon) lambda0 = -(lambda1 = 180);
  range[0] = lambda0, range[1] = lambda1;
  p0 = null;
} // Finds the left-right distance between two longitudes.
// This is almost the same as (lambda1 - lambda0 + 360°) % 360°, except that we want
// the distance between ±180° to be 360°.


function angle(lambda0, lambda1) {
  return (lambda1 -= lambda0) < 0 ? lambda1 + 360 : lambda1;
}

function rangeCompare(a, b) {
  return a[0] - b[0];
}

function rangeContains(range, x) {
  return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
}

function _default(feature) {
  var i, n, a, b, merged, deltaMax, delta;
  phi1 = lambda1 = -(lambda0 = phi0 = Infinity);
  ranges = [];
  (0, _stream.default)(feature, boundsStream); // First, sort ranges by their minimum longitudes.

  if (n = ranges.length) {
    ranges.sort(rangeCompare); // Then, merge any ranges that overlap.

    for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
      b = ranges[i];

      if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
        if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
        if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
      } else {
        merged.push(a = b);
      }
    } // Finally, find the largest gap between the merged ranges.
    // The final bounding box will be the inverse of this gap.


    for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
      b = merged[i];
      if ((delta = angle(a[1], b[0])) > deltaMax) deltaMax = delta, lambda0 = b[0], lambda1 = a[1];
    }
  }

  ranges = range = null;
  return lambda0 === Infinity || phi0 === Infinity ? [[NaN, NaN], [NaN, NaN]] : [[lambda0, phi0], [lambda1, phi1]];
}
},{"./adder":"../node_modules/d3-geo/src/adder.js","./area":"../node_modules/d3-geo/src/area.js","./cartesian":"../node_modules/d3-geo/src/cartesian.js","./math":"../node_modules/d3-geo/src/math.js","./stream":"../node_modules/d3-geo/src/stream.js"}],"../node_modules/d3-geo/src/centroid.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _math = require("./math");

var _noop = _interopRequireDefault(require("./noop"));

var _stream = _interopRequireDefault(require("./stream"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var W0, W1, X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, lambda00, phi00, // first point
x0, y0, z0; // previous point

var centroidStream = {
  sphere: _noop.default,
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function () {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function () {
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  }
}; // Arithmetic mean of Cartesian vectors.

function centroidPoint(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  var cosPhi = (0, _math.cos)(phi);
  centroidPointCartesian(cosPhi * (0, _math.cos)(lambda), cosPhi * (0, _math.sin)(lambda), (0, _math.sin)(phi));
}

function centroidPointCartesian(x, y, z) {
  ++W0;
  X0 += (x - X0) / W0;
  Y0 += (y - Y0) / W0;
  Z0 += (z - Z0) / W0;
}

function centroidLineStart() {
  centroidStream.point = centroidLinePointFirst;
}

function centroidLinePointFirst(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  var cosPhi = (0, _math.cos)(phi);
  x0 = cosPhi * (0, _math.cos)(lambda);
  y0 = cosPhi * (0, _math.sin)(lambda);
  z0 = (0, _math.sin)(phi);
  centroidStream.point = centroidLinePoint;
  centroidPointCartesian(x0, y0, z0);
}

function centroidLinePoint(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  var cosPhi = (0, _math.cos)(phi),
      x = cosPhi * (0, _math.cos)(lambda),
      y = cosPhi * (0, _math.sin)(lambda),
      z = (0, _math.sin)(phi),
      w = (0, _math.atan2)((0, _math.sqrt)((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w), x0 * x + y0 * y + z0 * z);
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

function centroidLineEnd() {
  centroidStream.point = centroidPoint;
} // See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
// J. Applied Mechanics 42, 239 (1975).


function centroidRingStart() {
  centroidStream.point = centroidRingPointFirst;
}

function centroidRingEnd() {
  centroidRingPoint(lambda00, phi00);
  centroidStream.point = centroidPoint;
}

function centroidRingPointFirst(lambda, phi) {
  lambda00 = lambda, phi00 = phi;
  lambda *= _math.radians, phi *= _math.radians;
  centroidStream.point = centroidRingPoint;
  var cosPhi = (0, _math.cos)(phi);
  x0 = cosPhi * (0, _math.cos)(lambda);
  y0 = cosPhi * (0, _math.sin)(lambda);
  z0 = (0, _math.sin)(phi);
  centroidPointCartesian(x0, y0, z0);
}

function centroidRingPoint(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  var cosPhi = (0, _math.cos)(phi),
      x = cosPhi * (0, _math.cos)(lambda),
      y = cosPhi * (0, _math.sin)(lambda),
      z = (0, _math.sin)(phi),
      cx = y0 * z - z0 * y,
      cy = z0 * x - x0 * z,
      cz = x0 * y - y0 * x,
      m = (0, _math.sqrt)(cx * cx + cy * cy + cz * cz),
      w = (0, _math.asin)(m),
      // line weight = angle
  v = m && -w / m; // area weight multiplier

  X2 += v * cx;
  Y2 += v * cy;
  Z2 += v * cz;
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

function _default(object) {
  W0 = W1 = X0 = Y0 = Z0 = X1 = Y1 = Z1 = X2 = Y2 = Z2 = 0;
  (0, _stream.default)(object, centroidStream);
  var x = X2,
      y = Y2,
      z = Z2,
      m = x * x + y * y + z * z; // If the area-weighted ccentroid is undefined, fall back to length-weighted ccentroid.

  if (m < _math.epsilon2) {
    x = X1, y = Y1, z = Z1; // If the feature has zero length, fall back to arithmetic mean of point vectors.

    if (W1 < _math.epsilon) x = X0, y = Y0, z = Z0;
    m = x * x + y * y + z * z; // If the feature still has an undefined ccentroid, then return.

    if (m < _math.epsilon2) return [NaN, NaN];
  }

  return [(0, _math.atan2)(y, x) * _math.degrees, (0, _math.asin)(z / (0, _math.sqrt)(m)) * _math.degrees];
}
},{"./math":"../node_modules/d3-geo/src/math.js","./noop":"../node_modules/d3-geo/src/noop.js","./stream":"../node_modules/d3-geo/src/stream.js"}],"../node_modules/d3-geo/src/constant.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return function () {
    return x;
  };
}
},{}],"../node_modules/d3-geo/src/compose.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(a, b) {
  function compose(x, y) {
    return x = a(x, y), b(x[0], x[1]);
  }

  if (a.invert && b.invert) compose.invert = function (x, y) {
    return x = b.invert(x, y), x && a.invert(x[0], x[1]);
  };
  return compose;
}
},{}],"../node_modules/d3-geo/src/rotation.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.rotateRadians = rotateRadians;
exports.default = _default;

var _compose = _interopRequireDefault(require("./compose"));

var _math = require("./math");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function rotationIdentity(lambda, phi) {
  return [(0, _math.abs)(lambda) > _math.pi ? lambda + Math.round(-lambda / _math.tau) * _math.tau : lambda, phi];
}

rotationIdentity.invert = rotationIdentity;

function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
  return (deltaLambda %= _math.tau) ? deltaPhi || deltaGamma ? (0, _compose.default)(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma)) : rotationLambda(deltaLambda) : deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma) : rotationIdentity;
}

function forwardRotationLambda(deltaLambda) {
  return function (lambda, phi) {
    return lambda += deltaLambda, [lambda > _math.pi ? lambda - _math.tau : lambda < -_math.pi ? lambda + _math.tau : lambda, phi];
  };
}

function rotationLambda(deltaLambda) {
  var rotation = forwardRotationLambda(deltaLambda);
  rotation.invert = forwardRotationLambda(-deltaLambda);
  return rotation;
}

function rotationPhiGamma(deltaPhi, deltaGamma) {
  var cosDeltaPhi = (0, _math.cos)(deltaPhi),
      sinDeltaPhi = (0, _math.sin)(deltaPhi),
      cosDeltaGamma = (0, _math.cos)(deltaGamma),
      sinDeltaGamma = (0, _math.sin)(deltaGamma);

  function rotation(lambda, phi) {
    var cosPhi = (0, _math.cos)(phi),
        x = (0, _math.cos)(lambda) * cosPhi,
        y = (0, _math.sin)(lambda) * cosPhi,
        z = (0, _math.sin)(phi),
        k = z * cosDeltaPhi + x * sinDeltaPhi;
    return [(0, _math.atan2)(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi), (0, _math.asin)(k * cosDeltaGamma + y * sinDeltaGamma)];
  }

  rotation.invert = function (lambda, phi) {
    var cosPhi = (0, _math.cos)(phi),
        x = (0, _math.cos)(lambda) * cosPhi,
        y = (0, _math.sin)(lambda) * cosPhi,
        z = (0, _math.sin)(phi),
        k = z * cosDeltaGamma - y * sinDeltaGamma;
    return [(0, _math.atan2)(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi), (0, _math.asin)(k * cosDeltaPhi - x * sinDeltaPhi)];
  };

  return rotation;
}

function _default(rotate) {
  rotate = rotateRadians(rotate[0] * _math.radians, rotate[1] * _math.radians, rotate.length > 2 ? rotate[2] * _math.radians : 0);

  function forward(coordinates) {
    coordinates = rotate(coordinates[0] * _math.radians, coordinates[1] * _math.radians);
    return coordinates[0] *= _math.degrees, coordinates[1] *= _math.degrees, coordinates;
  }

  forward.invert = function (coordinates) {
    coordinates = rotate.invert(coordinates[0] * _math.radians, coordinates[1] * _math.radians);
    return coordinates[0] *= _math.degrees, coordinates[1] *= _math.degrees, coordinates;
  };

  return forward;
}
},{"./compose":"../node_modules/d3-geo/src/compose.js","./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/circle.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.circleStream = circleStream;
exports.default = _default;

var _cartesian = require("./cartesian");

var _constant = _interopRequireDefault(require("./constant"));

var _math = require("./math");

var _rotation = require("./rotation");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

// Generates a circle centered at [0°, 0°], with a given radius and precision.
function circleStream(stream, radius, delta, direction, t0, t1) {
  if (!delta) return;
  var cosRadius = (0, _math.cos)(radius),
      sinRadius = (0, _math.sin)(radius),
      step = direction * delta;

  if (t0 == null) {
    t0 = radius + direction * _math.tau;
    t1 = radius - step / 2;
  } else {
    t0 = circleRadius(cosRadius, t0);
    t1 = circleRadius(cosRadius, t1);
    if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * _math.tau;
  }

  for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
    point = (0, _cartesian.spherical)([cosRadius, -sinRadius * (0, _math.cos)(t), -sinRadius * (0, _math.sin)(t)]);
    stream.point(point[0], point[1]);
  }
} // Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].


function circleRadius(cosRadius, point) {
  point = (0, _cartesian.cartesian)(point), point[0] -= cosRadius;
  (0, _cartesian.cartesianNormalizeInPlace)(point);
  var radius = (0, _math.acos)(-point[1]);
  return ((-point[2] < 0 ? -radius : radius) + _math.tau - _math.epsilon) % _math.tau;
}

function _default() {
  var center = (0, _constant.default)([0, 0]),
      radius = (0, _constant.default)(90),
      precision = (0, _constant.default)(6),
      ring,
      rotate,
      stream = {
    point: point
  };

  function point(x, y) {
    ring.push(x = rotate(x, y));
    x[0] *= _math.degrees, x[1] *= _math.degrees;
  }

  function circle() {
    var c = center.apply(this, arguments),
        r = radius.apply(this, arguments) * _math.radians,
        p = precision.apply(this, arguments) * _math.radians;

    ring = [];
    rotate = (0, _rotation.rotateRadians)(-c[0] * _math.radians, -c[1] * _math.radians, 0).invert;
    circleStream(stream, r, p, 1);
    c = {
      type: "Polygon",
      coordinates: [ring]
    };
    ring = rotate = null;
    return c;
  }

  circle.center = function (_) {
    return arguments.length ? (center = typeof _ === "function" ? _ : (0, _constant.default)([+_[0], +_[1]]), circle) : center;
  };

  circle.radius = function (_) {
    return arguments.length ? (radius = typeof _ === "function" ? _ : (0, _constant.default)(+_), circle) : radius;
  };

  circle.precision = function (_) {
    return arguments.length ? (precision = typeof _ === "function" ? _ : (0, _constant.default)(+_), circle) : precision;
  };

  return circle;
}
},{"./cartesian":"../node_modules/d3-geo/src/cartesian.js","./constant":"../node_modules/d3-geo/src/constant.js","./math":"../node_modules/d3-geo/src/math.js","./rotation":"../node_modules/d3-geo/src/rotation.js"}],"../node_modules/d3-geo/src/clip/buffer.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _noop = _interopRequireDefault(require("../noop"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  var lines = [],
      line;
  return {
    point: function (x, y) {
      line.push([x, y]);
    },
    lineStart: function () {
      lines.push(line = []);
    },
    lineEnd: _noop.default,
    rejoin: function () {
      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
    },
    result: function () {
      var result = lines;
      lines = [];
      line = null;
      return result;
    }
  };
}
},{"../noop":"../node_modules/d3-geo/src/noop.js"}],"../node_modules/d3-geo/src/pointEqual.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _math = require("./math");

function _default(a, b) {
  return (0, _math.abs)(a[0] - b[0]) < _math.epsilon && (0, _math.abs)(a[1] - b[1]) < _math.epsilon;
}
},{"./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/clip/rejoin.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _pointEqual = _interopRequireDefault(require("../pointEqual"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function Intersection(point, points, other, entry) {
  this.x = point;
  this.z = points;
  this.o = other; // another intersection

  this.e = entry; // is an entry?

  this.v = false; // visited

  this.n = this.p = null; // next & previous
} // A generalized polygon clipping algorithm: given a polygon that has been cut
// into its visible line segments, and rejoins the segments by interpolating
// along the clip edge.


function _default(segments, compareIntersection, startInside, interpolate, stream) {
  var subject = [],
      clip = [],
      i,
      n;
  segments.forEach(function (segment) {
    if ((n = segment.length - 1) <= 0) return;
    var n,
        p0 = segment[0],
        p1 = segment[n],
        x; // If the first and last points of a segment are coincident, then treat as a
    // closed ring. TODO if all rings are closed, then the winding order of the
    // exterior ring should be checked.

    if ((0, _pointEqual.default)(p0, p1)) {
      stream.lineStart();

      for (i = 0; i < n; ++i) stream.point((p0 = segment[i])[0], p0[1]);

      stream.lineEnd();
      return;
    }

    subject.push(x = new Intersection(p0, segment, null, true));
    clip.push(x.o = new Intersection(p0, null, x, false));
    subject.push(x = new Intersection(p1, segment, null, false));
    clip.push(x.o = new Intersection(p1, null, x, true));
  });
  if (!subject.length) return;
  clip.sort(compareIntersection);
  link(subject);
  link(clip);

  for (i = 0, n = clip.length; i < n; ++i) {
    clip[i].e = startInside = !startInside;
  }

  var start = subject[0],
      points,
      point;

  while (1) {
    // Find first unvisited intersection.
    var current = start,
        isSubject = true;

    while (current.v) if ((current = current.n) === start) return;

    points = current.z;
    stream.lineStart();

    do {
      current.v = current.o.v = true;

      if (current.e) {
        if (isSubject) {
          for (i = 0, n = points.length; i < n; ++i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.n.x, 1, stream);
        }

        current = current.n;
      } else {
        if (isSubject) {
          points = current.p.z;

          for (i = points.length - 1; i >= 0; --i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.p.x, -1, stream);
        }

        current = current.p;
      }

      current = current.o;
      points = current.z;
      isSubject = !isSubject;
    } while (!current.v);

    stream.lineEnd();
  }
}

function link(array) {
  if (!(n = array.length)) return;
  var n,
      i = 0,
      a = array[0],
      b;

  while (++i < n) {
    a.n = b = array[i];
    b.p = a;
    a = b;
  }

  a.n = b = array[0];
  b.p = a;
}
},{"../pointEqual":"../node_modules/d3-geo/src/pointEqual.js"}],"../node_modules/d3-geo/src/polygonContains.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _adder = _interopRequireDefault(require("./adder"));

var _cartesian = require("./cartesian");

var _math = require("./math");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var sum = (0, _adder.default)();

function _default(polygon, point) {
  var lambda = point[0],
      phi = point[1],
      sinPhi = (0, _math.sin)(phi),
      normal = [(0, _math.sin)(lambda), -(0, _math.cos)(lambda), 0],
      angle = 0,
      winding = 0;
  sum.reset();
  if (sinPhi === 1) phi = _math.halfPi + _math.epsilon;else if (sinPhi === -1) phi = -_math.halfPi - _math.epsilon;

  for (var i = 0, n = polygon.length; i < n; ++i) {
    if (!(m = (ring = polygon[i]).length)) continue;
    var ring,
        m,
        point0 = ring[m - 1],
        lambda0 = point0[0],
        phi0 = point0[1] / 2 + _math.quarterPi,
        sinPhi0 = (0, _math.sin)(phi0),
        cosPhi0 = (0, _math.cos)(phi0);

    for (var j = 0; j < m; ++j, lambda0 = lambda1, sinPhi0 = sinPhi1, cosPhi0 = cosPhi1, point0 = point1) {
      var point1 = ring[j],
          lambda1 = point1[0],
          phi1 = point1[1] / 2 + _math.quarterPi,
          sinPhi1 = (0, _math.sin)(phi1),
          cosPhi1 = (0, _math.cos)(phi1),
          delta = lambda1 - lambda0,
          sign = delta >= 0 ? 1 : -1,
          absDelta = sign * delta,
          antimeridian = absDelta > _math.pi,
          k = sinPhi0 * sinPhi1;
      sum.add((0, _math.atan2)(k * sign * (0, _math.sin)(absDelta), cosPhi0 * cosPhi1 + k * (0, _math.cos)(absDelta)));
      angle += antimeridian ? delta + sign * _math.tau : delta; // Are the longitudes either side of the point’s meridian (lambda),
      // and are the latitudes smaller than the parallel (phi)?

      if (antimeridian ^ lambda0 >= lambda ^ lambda1 >= lambda) {
        var arc = (0, _cartesian.cartesianCross)((0, _cartesian.cartesian)(point0), (0, _cartesian.cartesian)(point1));
        (0, _cartesian.cartesianNormalizeInPlace)(arc);
        var intersection = (0, _cartesian.cartesianCross)(normal, arc);
        (0, _cartesian.cartesianNormalizeInPlace)(intersection);
        var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * (0, _math.asin)(intersection[2]);

        if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
          winding += antimeridian ^ delta >= 0 ? 1 : -1;
        }
      }
    }
  } // First, determine whether the South pole is inside or outside:
  //
  // It is inside if:
  // * the polygon winds around it in a clockwise direction.
  // * the polygon does not (cumulatively) wind around it, but has a negative
  //   (counter-clockwise) area.
  //
  // Second, count the (signed) number of times a segment crosses a lambda
  // from the point to the South pole.  If it is zero, then the point is the
  // same side as the South pole.


  return (angle < -_math.epsilon || angle < _math.epsilon && sum < -_math.epsilon) ^ winding & 1;
}
},{"./adder":"../node_modules/d3-geo/src/adder.js","./cartesian":"../node_modules/d3-geo/src/cartesian.js","./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-array/src/ascending.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}
},{}],"../node_modules/d3-array/src/bisector.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _ascending = _interopRequireDefault(require("./ascending"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(compare) {
  if (compare.length === 1) compare = ascendingComparator(compare);
  return {
    left: function (a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;

      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) < 0) lo = mid + 1;else hi = mid;
      }

      return lo;
    },
    right: function (a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;

      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) > 0) hi = mid;else lo = mid + 1;
      }

      return lo;
    }
  };
}

function ascendingComparator(f) {
  return function (d, x) {
    return (0, _ascending.default)(f(d), x);
  };
}
},{"./ascending":"../node_modules/d3-array/src/ascending.js"}],"../node_modules/d3-array/src/bisect.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = exports.bisectLeft = exports.bisectRight = void 0;

var _ascending = _interopRequireDefault(require("./ascending"));

var _bisector = _interopRequireDefault(require("./bisector"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var ascendingBisect = (0, _bisector.default)(_ascending.default);
var bisectRight = ascendingBisect.right;
exports.bisectRight = bisectRight;
var bisectLeft = ascendingBisect.left;
exports.bisectLeft = bisectLeft;
var _default = bisectRight;
exports.default = _default;
},{"./ascending":"../node_modules/d3-array/src/ascending.js","./bisector":"../node_modules/d3-array/src/bisector.js"}],"../node_modules/d3-array/src/pairs.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.pair = pair;

function _default(array, f) {
  if (f == null) f = pair;
  var i = 0,
      n = array.length - 1,
      p = array[0],
      pairs = new Array(n < 0 ? 0 : n);

  while (i < n) pairs[i] = f(p, p = array[++i]);

  return pairs;
}

function pair(a, b) {
  return [a, b];
}
},{}],"../node_modules/d3-array/src/cross.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _pairs = require("./pairs");

function _default(values0, values1, reduce) {
  var n0 = values0.length,
      n1 = values1.length,
      values = new Array(n0 * n1),
      i0,
      i1,
      i,
      value0;
  if (reduce == null) reduce = _pairs.pair;

  for (i0 = i = 0; i0 < n0; ++i0) {
    for (value0 = values0[i0], i1 = 0; i1 < n1; ++i1, ++i) {
      values[i] = reduce(value0, values1[i1]);
    }
  }

  return values;
}
},{"./pairs":"../node_modules/d3-array/src/pairs.js"}],"../node_modules/d3-array/src/descending.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(a, b) {
  return b < a ? -1 : b > a ? 1 : b >= a ? 0 : NaN;
}
},{}],"../node_modules/d3-array/src/number.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return x === null ? NaN : +x;
}
},{}],"../node_modules/d3-array/src/variance.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _number = _interopRequireDefault(require("./number"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, valueof) {
  var n = values.length,
      m = 0,
      i = -1,
      mean = 0,
      value,
      delta,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(values[i]))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  } else {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(valueof(values[i], i, values)))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  if (m > 1) return sum / (m - 1);
}
},{"./number":"../node_modules/d3-array/src/number.js"}],"../node_modules/d3-array/src/deviation.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _variance = _interopRequireDefault(require("./variance"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(array, f) {
  var v = (0, _variance.default)(array, f);
  return v ? Math.sqrt(v) : v;
}
},{"./variance":"../node_modules/d3-array/src/variance.js"}],"../node_modules/d3-array/src/extent.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min,
      max;

  if (valueof == null) {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = max = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = values[i]) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  } else {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = max = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  return [min, max];
}
},{}],"../node_modules/d3-array/src/array.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.map = exports.slice = void 0;
var array = Array.prototype;
var slice = array.slice;
exports.slice = slice;
var map = array.map;
exports.map = map;
},{}],"../node_modules/d3-array/src/constant.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return function () {
    return x;
  };
}
},{}],"../node_modules/d3-array/src/identity.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return x;
}
},{}],"../node_modules/d3-array/src/range.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;
  var i = -1,
      n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
      range = new Array(n);

  while (++i < n) {
    range[i] = start + i * step;
  }

  return range;
}
},{}],"../node_modules/d3-array/src/ticks.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.tickIncrement = tickIncrement;
exports.tickStep = tickStep;
var e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

function _default(start, stop, count) {
  var reverse,
      i = -1,
      n,
      ticks,
      step;
  stop = +stop, start = +start, count = +count;
  if (start === stop && count > 0) return [start];
  if (reverse = stop < start) n = start, start = stop, stop = n;
  if ((step = tickIncrement(start, stop, count)) === 0 || !isFinite(step)) return [];

  if (step > 0) {
    start = Math.ceil(start / step);
    stop = Math.floor(stop / step);
    ticks = new Array(n = Math.ceil(stop - start + 1));

    while (++i < n) ticks[i] = (start + i) * step;
  } else {
    start = Math.floor(start * step);
    stop = Math.ceil(stop * step);
    ticks = new Array(n = Math.ceil(start - stop + 1));

    while (++i < n) ticks[i] = (start - i) / step;
  }

  if (reverse) ticks.reverse();
  return ticks;
}

function tickIncrement(start, stop, count) {
  var step = (stop - start) / Math.max(0, count),
      power = Math.floor(Math.log(step) / Math.LN10),
      error = step / Math.pow(10, power);
  return power >= 0 ? (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1) * Math.pow(10, power) : -Math.pow(10, -power) / (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1);
}

function tickStep(start, stop, count) {
  var step0 = Math.abs(stop - start) / Math.max(0, count),
      step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
      error = step0 / step1;
  if (error >= e10) step1 *= 10;else if (error >= e5) step1 *= 5;else if (error >= e2) step1 *= 2;
  return stop < start ? -step1 : step1;
}
},{}],"../node_modules/d3-array/src/threshold/sturges.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(values) {
  return Math.ceil(Math.log(values.length) / Math.LN2) + 1;
}
},{}],"../node_modules/d3-array/src/histogram.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _array = require("./array");

var _bisect = _interopRequireDefault(require("./bisect"));

var _constant = _interopRequireDefault(require("./constant"));

var _extent = _interopRequireDefault(require("./extent"));

var _identity = _interopRequireDefault(require("./identity"));

var _range = _interopRequireDefault(require("./range"));

var _ticks = require("./ticks");

var _sturges = _interopRequireDefault(require("./threshold/sturges"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  var value = _identity.default,
      domain = _extent.default,
      threshold = _sturges.default;

  function histogram(data) {
    var i,
        n = data.length,
        x,
        values = new Array(n);

    for (i = 0; i < n; ++i) {
      values[i] = value(data[i], i, data);
    }

    var xz = domain(values),
        x0 = xz[0],
        x1 = xz[1],
        tz = threshold(values, x0, x1); // Convert number of thresholds into uniform thresholds.

    if (!Array.isArray(tz)) {
      tz = (0, _ticks.tickStep)(x0, x1, tz);
      tz = (0, _range.default)(Math.ceil(x0 / tz) * tz, x1, tz); // exclusive
    } // Remove any thresholds outside the domain.


    var m = tz.length;

    while (tz[0] <= x0) tz.shift(), --m;

    while (tz[m - 1] > x1) tz.pop(), --m;

    var bins = new Array(m + 1),
        bin; // Initialize bins.

    for (i = 0; i <= m; ++i) {
      bin = bins[i] = [];
      bin.x0 = i > 0 ? tz[i - 1] : x0;
      bin.x1 = i < m ? tz[i] : x1;
    } // Assign data to bins by value, ignoring any outside the domain.


    for (i = 0; i < n; ++i) {
      x = values[i];

      if (x0 <= x && x <= x1) {
        bins[(0, _bisect.default)(tz, x, 0, m)].push(data[i]);
      }
    }

    return bins;
  }

  histogram.value = function (_) {
    return arguments.length ? (value = typeof _ === "function" ? _ : (0, _constant.default)(_), histogram) : value;
  };

  histogram.domain = function (_) {
    return arguments.length ? (domain = typeof _ === "function" ? _ : (0, _constant.default)([_[0], _[1]]), histogram) : domain;
  };

  histogram.thresholds = function (_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? (0, _constant.default)(_array.slice.call(_)) : (0, _constant.default)(_), histogram) : threshold;
  };

  return histogram;
}
},{"./array":"../node_modules/d3-array/src/array.js","./bisect":"../node_modules/d3-array/src/bisect.js","./constant":"../node_modules/d3-array/src/constant.js","./extent":"../node_modules/d3-array/src/extent.js","./identity":"../node_modules/d3-array/src/identity.js","./range":"../node_modules/d3-array/src/range.js","./ticks":"../node_modules/d3-array/src/ticks.js","./threshold/sturges":"../node_modules/d3-array/src/threshold/sturges.js"}],"../node_modules/d3-array/src/quantile.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _number = _interopRequireDefault(require("./number"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, p, valueof) {
  if (valueof == null) valueof = _number.default;
  if (!(n = values.length)) return;
  if ((p = +p) <= 0 || n < 2) return +valueof(values[0], 0, values);
  if (p >= 1) return +valueof(values[n - 1], n - 1, values);
  var n,
      i = (n - 1) * p,
      i0 = Math.floor(i),
      value0 = +valueof(values[i0], i0, values),
      value1 = +valueof(values[i0 + 1], i0 + 1, values);
  return value0 + (value1 - value0) * (i - i0);
}
},{"./number":"../node_modules/d3-array/src/number.js"}],"../node_modules/d3-array/src/threshold/freedmanDiaconis.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _array = require("../array");

var _ascending = _interopRequireDefault(require("../ascending"));

var _number = _interopRequireDefault(require("../number"));

var _quantile = _interopRequireDefault(require("../quantile"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, min, max) {
  values = _array.map.call(values, _number.default).sort(_ascending.default);
  return Math.ceil((max - min) / (2 * ((0, _quantile.default)(values, 0.75) - (0, _quantile.default)(values, 0.25)) * Math.pow(values.length, -1 / 3)));
}
},{"../array":"../node_modules/d3-array/src/array.js","../ascending":"../node_modules/d3-array/src/ascending.js","../number":"../node_modules/d3-array/src/number.js","../quantile":"../node_modules/d3-array/src/quantile.js"}],"../node_modules/d3-array/src/threshold/scott.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _deviation = _interopRequireDefault(require("../deviation"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, min, max) {
  return Math.ceil((max - min) / (3.5 * (0, _deviation.default)(values) * Math.pow(values.length, -1 / 3)));
}
},{"../deviation":"../node_modules/d3-array/src/deviation.js"}],"../node_modules/d3-array/src/max.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      max;

  if (valueof == null) {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        max = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = values[i]) != null && value > max) {
            max = value;
          }
        }
      }
    }
  } else {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        max = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  return max;
}
},{}],"../node_modules/d3-array/src/mean.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _number = _interopRequireDefault(require("./number"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, valueof) {
  var n = values.length,
      m = n,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(values[i]))) sum += value;else --m;
    }
  } else {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(valueof(values[i], i, values)))) sum += value;else --m;
    }
  }

  if (m) return sum / m;
}
},{"./number":"../node_modules/d3-array/src/number.js"}],"../node_modules/d3-array/src/median.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _ascending = _interopRequireDefault(require("./ascending"));

var _number = _interopRequireDefault(require("./number"));

var _quantile = _interopRequireDefault(require("./quantile"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      numbers = [];

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(values[i]))) {
        numbers.push(value);
      }
    }
  } else {
    while (++i < n) {
      if (!isNaN(value = (0, _number.default)(valueof(values[i], i, values)))) {
        numbers.push(value);
      }
    }
  }

  return (0, _quantile.default)(numbers.sort(_ascending.default), 0.5);
}
},{"./ascending":"../node_modules/d3-array/src/ascending.js","./number":"../node_modules/d3-array/src/number.js","./quantile":"../node_modules/d3-array/src/quantile.js"}],"../node_modules/d3-array/src/merge.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(arrays) {
  var n = arrays.length,
      m,
      i = -1,
      j = 0,
      merged,
      array;

  while (++i < n) j += arrays[i].length;

  merged = new Array(j);

  while (--n >= 0) {
    array = arrays[n];
    m = array.length;

    while (--m >= 0) {
      merged[--j] = array[m];
    }
  }

  return merged;
}
},{}],"../node_modules/d3-array/src/min.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min;

  if (valueof == null) {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = values[i]) != null && min > value) {
            min = value;
          }
        }
      }
    }
  } else {
    while (++i < n) {
      // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = value;

        while (++i < n) {
          // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  return min;
}
},{}],"../node_modules/d3-array/src/permute.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(array, indexes) {
  var i = indexes.length,
      permutes = new Array(i);

  while (i--) permutes[i] = array[indexes[i]];

  return permutes;
}
},{}],"../node_modules/d3-array/src/scan.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _ascending = _interopRequireDefault(require("./ascending"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(values, compare) {
  if (!(n = values.length)) return;
  var n,
      i = 0,
      j = 0,
      xi,
      xj = values[j];
  if (compare == null) compare = _ascending.default;

  while (++i < n) {
    if (compare(xi = values[i], xj) < 0 || compare(xj, xj) !== 0) {
      xj = xi, j = i;
    }
  }

  if (compare(xj, xj) === 0) return j;
}
},{"./ascending":"../node_modules/d3-array/src/ascending.js"}],"../node_modules/d3-array/src/shuffle.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(array, i0, i1) {
  var m = (i1 == null ? array.length : i1) - (i0 = i0 == null ? 0 : +i0),
      t,
      i;

  while (m) {
    i = Math.random() * m-- | 0;
    t = array[m + i0];
    array[m + i0] = array[i + i0];
    array[i + i0] = t;
  }

  return array;
}
},{}],"../node_modules/d3-array/src/sum.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (value = +values[i]) sum += value; // Note: zero and null are equivalent.
    }
  } else {
    while (++i < n) {
      if (value = +valueof(values[i], i, values)) sum += value;
    }
  }

  return sum;
}
},{}],"../node_modules/d3-array/src/transpose.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _min = _interopRequireDefault(require("./min"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(matrix) {
  if (!(n = matrix.length)) return [];

  for (var i = -1, m = (0, _min.default)(matrix, length), transpose = new Array(m); ++i < m;) {
    for (var j = -1, n, row = transpose[i] = new Array(n); ++j < n;) {
      row[j] = matrix[j][i];
    }
  }

  return transpose;
}

function length(d) {
  return d.length;
}
},{"./min":"../node_modules/d3-array/src/min.js"}],"../node_modules/d3-array/src/zip.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _transpose = _interopRequireDefault(require("./transpose"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  return (0, _transpose.default)(arguments);
}
},{"./transpose":"../node_modules/d3-array/src/transpose.js"}],"../node_modules/d3-array/src/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
Object.defineProperty(exports, "bisect", {
  enumerable: true,
  get: function () {
    return _bisect.default;
  }
});
Object.defineProperty(exports, "bisectRight", {
  enumerable: true,
  get: function () {
    return _bisect.bisectRight;
  }
});
Object.defineProperty(exports, "bisectLeft", {
  enumerable: true,
  get: function () {
    return _bisect.bisectLeft;
  }
});
Object.defineProperty(exports, "ascending", {
  enumerable: true,
  get: function () {
    return _ascending.default;
  }
});
Object.defineProperty(exports, "bisector", {
  enumerable: true,
  get: function () {
    return _bisector.default;
  }
});
Object.defineProperty(exports, "cross", {
  enumerable: true,
  get: function () {
    return _cross.default;
  }
});
Object.defineProperty(exports, "descending", {
  enumerable: true,
  get: function () {
    return _descending.default;
  }
});
Object.defineProperty(exports, "deviation", {
  enumerable: true,
  get: function () {
    return _deviation.default;
  }
});
Object.defineProperty(exports, "extent", {
  enumerable: true,
  get: function () {
    return _extent.default;
  }
});
Object.defineProperty(exports, "histogram", {
  enumerable: true,
  get: function () {
    return _histogram.default;
  }
});
Object.defineProperty(exports, "thresholdFreedmanDiaconis", {
  enumerable: true,
  get: function () {
    return _freedmanDiaconis.default;
  }
});
Object.defineProperty(exports, "thresholdScott", {
  enumerable: true,
  get: function () {
    return _scott.default;
  }
});
Object.defineProperty(exports, "thresholdSturges", {
  enumerable: true,
  get: function () {
    return _sturges.default;
  }
});
Object.defineProperty(exports, "max", {
  enumerable: true,
  get: function () {
    return _max.default;
  }
});
Object.defineProperty(exports, "mean", {
  enumerable: true,
  get: function () {
    return _mean.default;
  }
});
Object.defineProperty(exports, "median", {
  enumerable: true,
  get: function () {
    return _median.default;
  }
});
Object.defineProperty(exports, "merge", {
  enumerable: true,
  get: function () {
    return _merge.default;
  }
});
Object.defineProperty(exports, "min", {
  enumerable: true,
  get: function () {
    return _min.default;
  }
});
Object.defineProperty(exports, "pairs", {
  enumerable: true,
  get: function () {
    return _pairs.default;
  }
});
Object.defineProperty(exports, "permute", {
  enumerable: true,
  get: function () {
    return _permute.default;
  }
});
Object.defineProperty(exports, "quantile", {
  enumerable: true,
  get: function () {
    return _quantile.default;
  }
});
Object.defineProperty(exports, "range", {
  enumerable: true,
  get: function () {
    return _range.default;
  }
});
Object.defineProperty(exports, "scan", {
  enumerable: true,
  get: function () {
    return _scan.default;
  }
});
Object.defineProperty(exports, "shuffle", {
  enumerable: true,
  get: function () {
    return _shuffle.default;
  }
});
Object.defineProperty(exports, "sum", {
  enumerable: true,
  get: function () {
    return _sum.default;
  }
});
Object.defineProperty(exports, "ticks", {
  enumerable: true,
  get: function () {
    return _ticks.default;
  }
});
Object.defineProperty(exports, "tickIncrement", {
  enumerable: true,
  get: function () {
    return _ticks.tickIncrement;
  }
});
Object.defineProperty(exports, "tickStep", {
  enumerable: true,
  get: function () {
    return _ticks.tickStep;
  }
});
Object.defineProperty(exports, "transpose", {
  enumerable: true,
  get: function () {
    return _transpose.default;
  }
});
Object.defineProperty(exports, "variance", {
  enumerable: true,
  get: function () {
    return _variance.default;
  }
});
Object.defineProperty(exports, "zip", {
  enumerable: true,
  get: function () {
    return _zip.default;
  }
});

var _bisect = _interopRequireWildcard(require("./bisect"));

var _ascending = _interopRequireDefault(require("./ascending"));

var _bisector = _interopRequireDefault(require("./bisector"));

var _cross = _interopRequireDefault(require("./cross"));

var _descending = _interopRequireDefault(require("./descending"));

var _deviation = _interopRequireDefault(require("./deviation"));

var _extent = _interopRequireDefault(require("./extent"));

var _histogram = _interopRequireDefault(require("./histogram"));

var _freedmanDiaconis = _interopRequireDefault(require("./threshold/freedmanDiaconis"));

var _scott = _interopRequireDefault(require("./threshold/scott"));

var _sturges = _interopRequireDefault(require("./threshold/sturges"));

var _max = _interopRequireDefault(require("./max"));

var _mean = _interopRequireDefault(require("./mean"));

var _median = _interopRequireDefault(require("./median"));

var _merge = _interopRequireDefault(require("./merge"));

var _min = _interopRequireDefault(require("./min"));

var _pairs = _interopRequireDefault(require("./pairs"));

var _permute = _interopRequireDefault(require("./permute"));

var _quantile = _interopRequireDefault(require("./quantile"));

var _range = _interopRequireDefault(require("./range"));

var _scan = _interopRequireDefault(require("./scan"));

var _shuffle = _interopRequireDefault(require("./shuffle"));

var _sum = _interopRequireDefault(require("./sum"));

var _ticks = _interopRequireWildcard(require("./ticks"));

var _transpose = _interopRequireDefault(require("./transpose"));

var _variance = _interopRequireDefault(require("./variance"));

var _zip = _interopRequireDefault(require("./zip"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = Object.defineProperty && Object.getOwnPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : {}; if (desc.get || desc.set) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } } newObj.default = obj; return newObj; } }
},{"./bisect":"../node_modules/d3-array/src/bisect.js","./ascending":"../node_modules/d3-array/src/ascending.js","./bisector":"../node_modules/d3-array/src/bisector.js","./cross":"../node_modules/d3-array/src/cross.js","./descending":"../node_modules/d3-array/src/descending.js","./deviation":"../node_modules/d3-array/src/deviation.js","./extent":"../node_modules/d3-array/src/extent.js","./histogram":"../node_modules/d3-array/src/histogram.js","./threshold/freedmanDiaconis":"../node_modules/d3-array/src/threshold/freedmanDiaconis.js","./threshold/scott":"../node_modules/d3-array/src/threshold/scott.js","./threshold/sturges":"../node_modules/d3-array/src/threshold/sturges.js","./max":"../node_modules/d3-array/src/max.js","./mean":"../node_modules/d3-array/src/mean.js","./median":"../node_modules/d3-array/src/median.js","./merge":"../node_modules/d3-array/src/merge.js","./min":"../node_modules/d3-array/src/min.js","./pairs":"../node_modules/d3-array/src/pairs.js","./permute":"../node_modules/d3-array/src/permute.js","./quantile":"../node_modules/d3-array/src/quantile.js","./range":"../node_modules/d3-array/src/range.js","./scan":"../node_modules/d3-array/src/scan.js","./shuffle":"../node_modules/d3-array/src/shuffle.js","./sum":"../node_modules/d3-array/src/sum.js","./ticks":"../node_modules/d3-array/src/ticks.js","./transpose":"../node_modules/d3-array/src/transpose.js","./variance":"../node_modules/d3-array/src/variance.js","./zip":"../node_modules/d3-array/src/zip.js"}],"../node_modules/d3-geo/src/clip/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _buffer = _interopRequireDefault(require("./buffer"));

var _rejoin = _interopRequireDefault(require("./rejoin"));

var _math = require("../math");

var _polygonContains = _interopRequireDefault(require("../polygonContains"));

var _d3Array = require("d3-array");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(pointVisible, clipLine, interpolate, start) {
  return function (sink) {
    var line = clipLine(sink),
        ringBuffer = (0, _buffer.default)(),
        ringSink = clipLine(ringBuffer),
        polygonStarted = false,
        polygon,
        segments,
        ring;
    var clip = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function () {
        clip.point = pointRing;
        clip.lineStart = ringStart;
        clip.lineEnd = ringEnd;
        segments = [];
        polygon = [];
      },
      polygonEnd: function () {
        clip.point = point;
        clip.lineStart = lineStart;
        clip.lineEnd = lineEnd;
        segments = (0, _d3Array.merge)(segments);
        var startInside = (0, _polygonContains.default)(polygon, start);

        if (segments.length) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          (0, _rejoin.default)(segments, compareIntersection, startInside, interpolate, sink);
        } else if (startInside) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
        }

        if (polygonStarted) sink.polygonEnd(), polygonStarted = false;
        segments = polygon = null;
      },
      sphere: function () {
        sink.polygonStart();
        sink.lineStart();
        interpolate(null, null, 1, sink);
        sink.lineEnd();
        sink.polygonEnd();
      }
    };

    function point(lambda, phi) {
      if (pointVisible(lambda, phi)) sink.point(lambda, phi);
    }

    function pointLine(lambda, phi) {
      line.point(lambda, phi);
    }

    function lineStart() {
      clip.point = pointLine;
      line.lineStart();
    }

    function lineEnd() {
      clip.point = point;
      line.lineEnd();
    }

    function pointRing(lambda, phi) {
      ring.push([lambda, phi]);
      ringSink.point(lambda, phi);
    }

    function ringStart() {
      ringSink.lineStart();
      ring = [];
    }

    function ringEnd() {
      pointRing(ring[0][0], ring[0][1]);
      ringSink.lineEnd();
      var clean = ringSink.clean(),
          ringSegments = ringBuffer.result(),
          i,
          n = ringSegments.length,
          m,
          segment,
          point;
      ring.pop();
      polygon.push(ring);
      ring = null;
      if (!n) return; // No intersections.

      if (clean & 1) {
        segment = ringSegments[0];

        if ((m = segment.length - 1) > 0) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();

          for (i = 0; i < m; ++i) sink.point((point = segment[i])[0], point[1]);

          sink.lineEnd();
        }

        return;
      } // Rejoin connected segments.
      // TODO reuse ringBuffer.rejoin()?


      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));
      segments.push(ringSegments.filter(validSegment));
    }

    return clip;
  };
}

function validSegment(segment) {
  return segment.length > 1;
} // Intersections are sorted along the clip edge. For both antimeridian cutting
// and circle clipping, the same comparison is used.


function compareIntersection(a, b) {
  return ((a = a.x)[0] < 0 ? a[1] - _math.halfPi - _math.epsilon : _math.halfPi - a[1]) - ((b = b.x)[0] < 0 ? b[1] - _math.halfPi - _math.epsilon : _math.halfPi - b[1]);
}
},{"./buffer":"../node_modules/d3-geo/src/clip/buffer.js","./rejoin":"../node_modules/d3-geo/src/clip/rejoin.js","../math":"../node_modules/d3-geo/src/math.js","../polygonContains":"../node_modules/d3-geo/src/polygonContains.js","d3-array":"../node_modules/d3-array/src/index.js"}],"../node_modules/d3-geo/src/clip/antimeridian.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _index = _interopRequireDefault(require("./index"));

var _math = require("../math");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var _default = (0, _index.default)(function () {
  return true;
}, clipAntimeridianLine, clipAntimeridianInterpolate, [-_math.pi, -_math.halfPi]); // Takes a line and cuts into visible segments. Return values: 0 - there were
// intersections or the line was empty; 1 - no intersections; 2 - there were
// intersections, and the first and last segments should be rejoined.


exports.default = _default;

function clipAntimeridianLine(stream) {
  var lambda0 = NaN,
      phi0 = NaN,
      sign0 = NaN,
      clean; // no intersections

  return {
    lineStart: function () {
      stream.lineStart();
      clean = 1;
    },
    point: function (lambda1, phi1) {
      var sign1 = lambda1 > 0 ? _math.pi : -_math.pi,
          delta = (0, _math.abs)(lambda1 - lambda0);

      if ((0, _math.abs)(delta - _math.pi) < _math.epsilon) {
        // line crosses a pole
        stream.point(lambda0, phi0 = (phi0 + phi1) / 2 > 0 ? _math.halfPi : -_math.halfPi);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        stream.point(lambda1, phi0);
        clean = 0;
      } else if (sign0 !== sign1 && delta >= _math.pi) {
        // line crosses antimeridian
        if ((0, _math.abs)(lambda0 - sign0) < _math.epsilon) lambda0 -= sign0 * _math.epsilon; // handle degeneracies

        if ((0, _math.abs)(lambda1 - sign1) < _math.epsilon) lambda1 -= sign1 * _math.epsilon;
        phi0 = clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        clean = 0;
      }

      stream.point(lambda0 = lambda1, phi0 = phi1);
      sign0 = sign1;
    },
    lineEnd: function () {
      stream.lineEnd();
      lambda0 = phi0 = NaN;
    },
    clean: function () {
      return 2 - clean; // if intersections, rejoin first and last segments
    }
  };
}

function clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1) {
  var cosPhi0,
      cosPhi1,
      sinLambda0Lambda1 = (0, _math.sin)(lambda0 - lambda1);
  return (0, _math.abs)(sinLambda0Lambda1) > _math.epsilon ? (0, _math.atan)(((0, _math.sin)(phi0) * (cosPhi1 = (0, _math.cos)(phi1)) * (0, _math.sin)(lambda1) - (0, _math.sin)(phi1) * (cosPhi0 = (0, _math.cos)(phi0)) * (0, _math.sin)(lambda0)) / (cosPhi0 * cosPhi1 * sinLambda0Lambda1)) : (phi0 + phi1) / 2;
}

function clipAntimeridianInterpolate(from, to, direction, stream) {
  var phi;

  if (from == null) {
    phi = direction * _math.halfPi;
    stream.point(-_math.pi, phi);
    stream.point(0, phi);
    stream.point(_math.pi, phi);
    stream.point(_math.pi, 0);
    stream.point(_math.pi, -phi);
    stream.point(0, -phi);
    stream.point(-_math.pi, -phi);
    stream.point(-_math.pi, 0);
    stream.point(-_math.pi, phi);
  } else if ((0, _math.abs)(from[0] - to[0]) > _math.epsilon) {
    var lambda = from[0] < to[0] ? _math.pi : -_math.pi;
    phi = direction * lambda / 2;
    stream.point(-lambda, phi);
    stream.point(0, phi);
    stream.point(lambda, phi);
  } else {
    stream.point(to[0], to[1]);
  }
}
},{"./index":"../node_modules/d3-geo/src/clip/index.js","../math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/clip/circle.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _cartesian = require("../cartesian");

var _circle = require("../circle");

var _math = require("../math");

var _pointEqual = _interopRequireDefault(require("../pointEqual"));

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(radius) {
  var cr = (0, _math.cos)(radius),
      delta = 6 * _math.radians,
      smallRadius = cr > 0,
      notHemisphere = (0, _math.abs)(cr) > _math.epsilon; // TODO optimise for this common case


  function interpolate(from, to, direction, stream) {
    (0, _circle.circleStream)(stream, radius, delta, direction, from, to);
  }

  function visible(lambda, phi) {
    return (0, _math.cos)(lambda) * (0, _math.cos)(phi) > cr;
  } // Takes a line and cuts into visible segments. Return values used for polygon
  // clipping: 0 - there were intersections or the line was empty; 1 - no
  // intersections 2 - there were intersections, and the first and last segments
  // should be rejoined.


  function clipLine(stream) {
    var point0, // previous point
    c0, // code for previous point
    v0, // visibility of previous point
    v00, // visibility of first point
    clean; // no intersections

    return {
      lineStart: function () {
        v00 = v0 = false;
        clean = 1;
      },
      point: function (lambda, phi) {
        var point1 = [lambda, phi],
            point2,
            v = visible(lambda, phi),
            c = smallRadius ? v ? 0 : code(lambda, phi) : v ? code(lambda + (lambda < 0 ? _math.pi : -_math.pi), phi) : 0;
        if (!point0 && (v00 = v0 = v)) stream.lineStart(); // Handle degeneracies.
        // TODO ignore if not clipping polygons.

        if (v !== v0) {
          point2 = intersect(point0, point1);

          if (!point2 || (0, _pointEqual.default)(point0, point2) || (0, _pointEqual.default)(point1, point2)) {
            point1[0] += _math.epsilon;
            point1[1] += _math.epsilon;
            v = visible(point1[0], point1[1]);
          }
        }

        if (v !== v0) {
          clean = 0;

          if (v) {
            // outside going in
            stream.lineStart();
            point2 = intersect(point1, point0);
            stream.point(point2[0], point2[1]);
          } else {
            // inside going out
            point2 = intersect(point0, point1);
            stream.point(point2[0], point2[1]);
            stream.lineEnd();
          }

          point0 = point2;
        } else if (notHemisphere && point0 && smallRadius ^ v) {
          var t; // If the codes for two points are different, or are both zero,
          // and there this segment intersects with the small circle.

          if (!(c & c0) && (t = intersect(point1, point0, true))) {
            clean = 0;

            if (smallRadius) {
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
            } else {
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
            }
          }
        }

        if (v && (!point0 || !(0, _pointEqual.default)(point0, point1))) {
          stream.point(point1[0], point1[1]);
        }

        point0 = point1, v0 = v, c0 = c;
      },
      lineEnd: function () {
        if (v0) stream.lineEnd();
        point0 = null;
      },
      // Rejoin first and last segments if there were intersections and the first
      // and last points were visible.
      clean: function () {
        return clean | (v00 && v0) << 1;
      }
    };
  } // Intersects the great circle between a and b with the clip circle.


  function intersect(a, b, two) {
    var pa = (0, _cartesian.cartesian)(a),
        pb = (0, _cartesian.cartesian)(b); // We have two planes, n1.p = d1 and n2.p = d2.
    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).

    var n1 = [1, 0, 0],
        // normal
    n2 = (0, _cartesian.cartesianCross)(pa, pb),
        n2n2 = (0, _cartesian.cartesianDot)(n2, n2),
        n1n2 = n2[0],
        // cartesianDot(n1, n2),
    determinant = n2n2 - n1n2 * n1n2; // Two polar points.

    if (!determinant) return !two && a;
    var c1 = cr * n2n2 / determinant,
        c2 = -cr * n1n2 / determinant,
        n1xn2 = (0, _cartesian.cartesianCross)(n1, n2),
        A = (0, _cartesian.cartesianScale)(n1, c1),
        B = (0, _cartesian.cartesianScale)(n2, c2);
    (0, _cartesian.cartesianAddInPlace)(A, B); // Solve |p(t)|^2 = 1.

    var u = n1xn2,
        w = (0, _cartesian.cartesianDot)(A, u),
        uu = (0, _cartesian.cartesianDot)(u, u),
        t2 = w * w - uu * ((0, _cartesian.cartesianDot)(A, A) - 1);
    if (t2 < 0) return;
    var t = (0, _math.sqrt)(t2),
        q = (0, _cartesian.cartesianScale)(u, (-w - t) / uu);
    (0, _cartesian.cartesianAddInPlace)(q, A);
    q = (0, _cartesian.spherical)(q);
    if (!two) return q; // Two intersection points.

    var lambda0 = a[0],
        lambda1 = b[0],
        phi0 = a[1],
        phi1 = b[1],
        z;
    if (lambda1 < lambda0) z = lambda0, lambda0 = lambda1, lambda1 = z;

    var delta = lambda1 - lambda0,
        polar = (0, _math.abs)(delta - _math.pi) < _math.epsilon,
        meridian = polar || delta < _math.epsilon;

    if (!polar && phi1 < phi0) z = phi0, phi0 = phi1, phi1 = z; // Check that the first point is between a and b.

    if (meridian ? polar ? phi0 + phi1 > 0 ^ q[1] < ((0, _math.abs)(q[0] - lambda0) < _math.epsilon ? phi0 : phi1) : phi0 <= q[1] && q[1] <= phi1 : delta > _math.pi ^ (lambda0 <= q[0] && q[0] <= lambda1)) {
      var q1 = (0, _cartesian.cartesianScale)(u, (-w + t) / uu);
      (0, _cartesian.cartesianAddInPlace)(q1, A);
      return [q, (0, _cartesian.spherical)(q1)];
    }
  } // Generates a 4-bit vector representing the location of a point relative to
  // the small circle's bounding box.


  function code(lambda, phi) {
    var r = smallRadius ? radius : _math.pi - radius,
        code = 0;
    if (lambda < -r) code |= 1; // left
    else if (lambda > r) code |= 2; // right

    if (phi < -r) code |= 4; // below
    else if (phi > r) code |= 8; // above

    return code;
  }

  return (0, _index.default)(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-_math.pi, radius - _math.pi]);
}
},{"../cartesian":"../node_modules/d3-geo/src/cartesian.js","../circle":"../node_modules/d3-geo/src/circle.js","../math":"../node_modules/d3-geo/src/math.js","../pointEqual":"../node_modules/d3-geo/src/pointEqual.js","./index":"../node_modules/d3-geo/src/clip/index.js"}],"../node_modules/d3-geo/src/clip/line.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(a, b, x0, y0, x1, y1) {
  var ax = a[0],
      ay = a[1],
      bx = b[0],
      by = b[1],
      t0 = 0,
      t1 = 1,
      dx = bx - ax,
      dy = by - ay,
      r;
  r = x0 - ax;
  if (!dx && r > 0) return;
  r /= dx;

  if (dx < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dx > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = x1 - ax;
  if (!dx && r < 0) return;
  r /= dx;

  if (dx < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dx > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  r = y0 - ay;
  if (!dy && r > 0) return;
  r /= dy;

  if (dy < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dy > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = y1 - ay;
  if (!dy && r < 0) return;
  r /= dy;

  if (dy < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dy > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
  if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
  return true;
}
},{}],"../node_modules/d3-geo/src/clip/rectangle.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = clipRectangle;

var _math = require("../math");

var _buffer = _interopRequireDefault(require("./buffer"));

var _line = _interopRequireDefault(require("./line"));

var _rejoin = _interopRequireDefault(require("./rejoin"));

var _d3Array = require("d3-array");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var clipMax = 1e9,
    clipMin = -clipMax; // TODO Use d3-polygon’s polygonContains here for the ring check?
// TODO Eliminate duplicate buffering in clipBuffer and polygon.push?

function clipRectangle(x0, y0, x1, y1) {
  function visible(x, y) {
    return x0 <= x && x <= x1 && y0 <= y && y <= y1;
  }

  function interpolate(from, to, direction, stream) {
    var a = 0,
        a1 = 0;

    if (from == null || (a = corner(from, direction)) !== (a1 = corner(to, direction)) || comparePoint(from, to) < 0 ^ direction > 0) {
      do stream.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0); while ((a = (a + direction + 4) % 4) !== a1);
    } else {
      stream.point(to[0], to[1]);
    }
  }

  function corner(p, direction) {
    return (0, _math.abs)(p[0] - x0) < _math.epsilon ? direction > 0 ? 0 : 3 : (0, _math.abs)(p[0] - x1) < _math.epsilon ? direction > 0 ? 2 : 1 : (0, _math.abs)(p[1] - y0) < _math.epsilon ? direction > 0 ? 1 : 0 : direction > 0 ? 3 : 2; // abs(p[1] - y1) < epsilon
  }

  function compareIntersection(a, b) {
    return comparePoint(a.x, b.x);
  }

  function comparePoint(a, b) {
    var ca = corner(a, 1),
        cb = corner(b, 1);
    return ca !== cb ? ca - cb : ca === 0 ? b[1] - a[1] : ca === 1 ? a[0] - b[0] : ca === 2 ? a[1] - b[1] : b[0] - a[0];
  }

  return function (stream) {
    var activeStream = stream,
        bufferStream = (0, _buffer.default)(),
        segments,
        polygon,
        ring,
        x__,
        y__,
        v__,
        // first point
    x_,
        y_,
        v_,
        // previous point
    first,
        clean;
    var clipStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: polygonStart,
      polygonEnd: polygonEnd
    };

    function point(x, y) {
      if (visible(x, y)) activeStream.point(x, y);
    }

    function polygonInside() {
      var winding = 0;

      for (var i = 0, n = polygon.length; i < n; ++i) {
        for (var ring = polygon[i], j = 1, m = ring.length, point = ring[0], a0, a1, b0 = point[0], b1 = point[1]; j < m; ++j) {
          a0 = b0, a1 = b1, point = ring[j], b0 = point[0], b1 = point[1];

          if (a1 <= y1) {
            if (b1 > y1 && (b0 - a0) * (y1 - a1) > (b1 - a1) * (x0 - a0)) ++winding;
          } else {
            if (b1 <= y1 && (b0 - a0) * (y1 - a1) < (b1 - a1) * (x0 - a0)) --winding;
          }
        }
      }

      return winding;
    } // Buffer geometry within a polygon and then clip it en masse.


    function polygonStart() {
      activeStream = bufferStream, segments = [], polygon = [], clean = true;
    }

    function polygonEnd() {
      var startInside = polygonInside(),
          cleanInside = clean && startInside,
          visible = (segments = (0, _d3Array.merge)(segments)).length;

      if (cleanInside || visible) {
        stream.polygonStart();

        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }

        if (visible) {
          (0, _rejoin.default)(segments, compareIntersection, startInside, interpolate, stream);
        }

        stream.polygonEnd();
      }

      activeStream = stream, segments = polygon = ring = null;
    }

    function lineStart() {
      clipStream.point = linePoint;
      if (polygon) polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = NaN;
    } // TODO rather than special-case polygons, simply handle them separately.
    // Ideally, coincident intersection points should be jittered to avoid
    // clipping issues.


    function lineEnd() {
      if (segments) {
        linePoint(x__, y__);
        if (v__ && v_) bufferStream.rejoin();
        segments.push(bufferStream.result());
      }

      clipStream.point = point;
      if (v_) activeStream.lineEnd();
    }

    function linePoint(x, y) {
      var v = visible(x, y);
      if (polygon) ring.push([x, y]);

      if (first) {
        x__ = x, y__ = y, v__ = v;
        first = false;

        if (v) {
          activeStream.lineStart();
          activeStream.point(x, y);
        }
      } else {
        if (v && v_) activeStream.point(x, y);else {
          var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))],
              b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];

          if ((0, _line.default)(a, b, x0, y0, x1, y1)) {
            if (!v_) {
              activeStream.lineStart();
              activeStream.point(a[0], a[1]);
            }

            activeStream.point(b[0], b[1]);
            if (!v) activeStream.lineEnd();
            clean = false;
          } else if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
            clean = false;
          }
        }
      }

      x_ = x, y_ = y, v_ = v;
    }

    return clipStream;
  };
}
},{"../math":"../node_modules/d3-geo/src/math.js","./buffer":"../node_modules/d3-geo/src/clip/buffer.js","./line":"../node_modules/d3-geo/src/clip/line.js","./rejoin":"../node_modules/d3-geo/src/clip/rejoin.js","d3-array":"../node_modules/d3-array/src/index.js"}],"../node_modules/d3-geo/src/clip/extent.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _rectangle = _interopRequireDefault(require("./rectangle"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  var x0 = 0,
      y0 = 0,
      x1 = 960,
      y1 = 500,
      cache,
      cacheStream,
      clip;
  return clip = {
    stream: function (stream) {
      return cache && cacheStream === stream ? cache : cache = (0, _rectangle.default)(x0, y0, x1, y1)(cacheStream = stream);
    },
    extent: function (_) {
      return arguments.length ? (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1], cache = cacheStream = null, clip) : [[x0, y0], [x1, y1]];
    }
  };
}
},{"./rectangle":"../node_modules/d3-geo/src/clip/rectangle.js"}],"../node_modules/d3-geo/src/length.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _adder = _interopRequireDefault(require("./adder"));

var _math = require("./math");

var _noop = _interopRequireDefault(require("./noop"));

var _stream = _interopRequireDefault(require("./stream"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var lengthSum = (0, _adder.default)(),
    lambda0,
    sinPhi0,
    cosPhi0;
var lengthStream = {
  sphere: _noop.default,
  point: _noop.default,
  lineStart: lengthLineStart,
  lineEnd: _noop.default,
  polygonStart: _noop.default,
  polygonEnd: _noop.default
};

function lengthLineStart() {
  lengthStream.point = lengthPointFirst;
  lengthStream.lineEnd = lengthLineEnd;
}

function lengthLineEnd() {
  lengthStream.point = lengthStream.lineEnd = _noop.default;
}

function lengthPointFirst(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  lambda0 = lambda, sinPhi0 = (0, _math.sin)(phi), cosPhi0 = (0, _math.cos)(phi);
  lengthStream.point = lengthPoint;
}

function lengthPoint(lambda, phi) {
  lambda *= _math.radians, phi *= _math.radians;
  var sinPhi = (0, _math.sin)(phi),
      cosPhi = (0, _math.cos)(phi),
      delta = (0, _math.abs)(lambda - lambda0),
      cosDelta = (0, _math.cos)(delta),
      sinDelta = (0, _math.sin)(delta),
      x = cosPhi * sinDelta,
      y = cosPhi0 * sinPhi - sinPhi0 * cosPhi * cosDelta,
      z = sinPhi0 * sinPhi + cosPhi0 * cosPhi * cosDelta;
  lengthSum.add((0, _math.atan2)((0, _math.sqrt)(x * x + y * y), z));
  lambda0 = lambda, sinPhi0 = sinPhi, cosPhi0 = cosPhi;
}

function _default(object) {
  lengthSum.reset();
  (0, _stream.default)(object, lengthStream);
  return +lengthSum;
}
},{"./adder":"../node_modules/d3-geo/src/adder.js","./math":"../node_modules/d3-geo/src/math.js","./noop":"../node_modules/d3-geo/src/noop.js","./stream":"../node_modules/d3-geo/src/stream.js"}],"../node_modules/d3-geo/src/distance.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _length = _interopRequireDefault(require("./length"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var coordinates = [null, null],
    object = {
  type: "LineString",
  coordinates: coordinates
};

function _default(a, b) {
  coordinates[0] = a;
  coordinates[1] = b;
  return (0, _length.default)(object);
}
},{"./length":"../node_modules/d3-geo/src/length.js"}],"../node_modules/d3-geo/src/contains.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _polygonContains = _interopRequireDefault(require("./polygonContains"));

var _distance = _interopRequireDefault(require("./distance"));

var _math = require("./math");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var containsObjectType = {
  Feature: function (object, point) {
    return containsGeometry(object.geometry, point);
  },
  FeatureCollection: function (object, point) {
    var features = object.features,
        i = -1,
        n = features.length;

    while (++i < n) if (containsGeometry(features[i].geometry, point)) return true;

    return false;
  }
};
var containsGeometryType = {
  Sphere: function () {
    return true;
  },
  Point: function (object, point) {
    return containsPoint(object.coordinates, point);
  },
  MultiPoint: function (object, point) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) if (containsPoint(coordinates[i], point)) return true;

    return false;
  },
  LineString: function (object, point) {
    return containsLine(object.coordinates, point);
  },
  MultiLineString: function (object, point) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) if (containsLine(coordinates[i], point)) return true;

    return false;
  },
  Polygon: function (object, point) {
    return containsPolygon(object.coordinates, point);
  },
  MultiPolygon: function (object, point) {
    var coordinates = object.coordinates,
        i = -1,
        n = coordinates.length;

    while (++i < n) if (containsPolygon(coordinates[i], point)) return true;

    return false;
  },
  GeometryCollection: function (object, point) {
    var geometries = object.geometries,
        i = -1,
        n = geometries.length;

    while (++i < n) if (containsGeometry(geometries[i], point)) return true;

    return false;
  }
};

function containsGeometry(geometry, point) {
  return geometry && containsGeometryType.hasOwnProperty(geometry.type) ? containsGeometryType[geometry.type](geometry, point) : false;
}

function containsPoint(coordinates, point) {
  return (0, _distance.default)(coordinates, point) === 0;
}

function containsLine(coordinates, point) {
  var ab = (0, _distance.default)(coordinates[0], coordinates[1]),
      ao = (0, _distance.default)(coordinates[0], point),
      ob = (0, _distance.default)(point, coordinates[1]);
  return ao + ob <= ab + _math.epsilon;
}

function containsPolygon(coordinates, point) {
  return !!(0, _polygonContains.default)(coordinates.map(ringRadians), pointRadians(point));
}

function ringRadians(ring) {
  return ring = ring.map(pointRadians), ring.pop(), ring;
}

function pointRadians(point) {
  return [point[0] * _math.radians, point[1] * _math.radians];
}

function _default(object, point) {
  return (object && containsObjectType.hasOwnProperty(object.type) ? containsObjectType[object.type] : containsGeometry)(object, point);
}
},{"./polygonContains":"../node_modules/d3-geo/src/polygonContains.js","./distance":"../node_modules/d3-geo/src/distance.js","./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/graticule.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = graticule;
exports.graticule10 = graticule10;

var _d3Array = require("d3-array");

var _math = require("./math");

function graticuleX(y0, y1, dy) {
  var y = (0, _d3Array.range)(y0, y1 - _math.epsilon, dy).concat(y1);
  return function (x) {
    return y.map(function (y) {
      return [x, y];
    });
  };
}

function graticuleY(x0, x1, dx) {
  var x = (0, _d3Array.range)(x0, x1 - _math.epsilon, dx).concat(x1);
  return function (y) {
    return x.map(function (x) {
      return [x, y];
    });
  };
}

function graticule() {
  var x1,
      x0,
      X1,
      X0,
      y1,
      y0,
      Y1,
      Y0,
      dx = 10,
      dy = dx,
      DX = 90,
      DY = 360,
      x,
      y,
      X,
      Y,
      precision = 2.5;

  function graticule() {
    return {
      type: "MultiLineString",
      coordinates: lines()
    };
  }

  function lines() {
    return (0, _d3Array.range)((0, _math.ceil)(X0 / DX) * DX, X1, DX).map(X).concat((0, _d3Array.range)((0, _math.ceil)(Y0 / DY) * DY, Y1, DY).map(Y)).concat((0, _d3Array.range)((0, _math.ceil)(x0 / dx) * dx, x1, dx).filter(function (x) {
      return (0, _math.abs)(x % DX) > _math.epsilon;
    }).map(x)).concat((0, _d3Array.range)((0, _math.ceil)(y0 / dy) * dy, y1, dy).filter(function (y) {
      return (0, _math.abs)(y % DY) > _math.epsilon;
    }).map(y));
  }

  graticule.lines = function () {
    return lines().map(function (coordinates) {
      return {
        type: "LineString",
        coordinates: coordinates
      };
    });
  };

  graticule.outline = function () {
    return {
      type: "Polygon",
      coordinates: [X(X0).concat(Y(Y1).slice(1), X(X1).reverse().slice(1), Y(Y0).reverse().slice(1))]
    };
  };

  graticule.extent = function (_) {
    if (!arguments.length) return graticule.extentMinor();
    return graticule.extentMajor(_).extentMinor(_);
  };

  graticule.extentMajor = function (_) {
    if (!arguments.length) return [[X0, Y0], [X1, Y1]];
    X0 = +_[0][0], X1 = +_[1][0];
    Y0 = +_[0][1], Y1 = +_[1][1];
    if (X0 > X1) _ = X0, X0 = X1, X1 = _;
    if (Y0 > Y1) _ = Y0, Y0 = Y1, Y1 = _;
    return graticule.precision(precision);
  };

  graticule.extentMinor = function (_) {
    if (!arguments.length) return [[x0, y0], [x1, y1]];
    x0 = +_[0][0], x1 = +_[1][0];
    y0 = +_[0][1], y1 = +_[1][1];
    if (x0 > x1) _ = x0, x0 = x1, x1 = _;
    if (y0 > y1) _ = y0, y0 = y1, y1 = _;
    return graticule.precision(precision);
  };

  graticule.step = function (_) {
    if (!arguments.length) return graticule.stepMinor();
    return graticule.stepMajor(_).stepMinor(_);
  };

  graticule.stepMajor = function (_) {
    if (!arguments.length) return [DX, DY];
    DX = +_[0], DY = +_[1];
    return graticule;
  };

  graticule.stepMinor = function (_) {
    if (!arguments.length) return [dx, dy];
    dx = +_[0], dy = +_[1];
    return graticule;
  };

  graticule.precision = function (_) {
    if (!arguments.length) return precision;
    precision = +_;
    x = graticuleX(y0, y1, 90);
    y = graticuleY(x0, x1, precision);
    X = graticuleX(Y0, Y1, 90);
    Y = graticuleY(X0, X1, precision);
    return graticule;
  };

  return graticule.extentMajor([[-180, -90 + _math.epsilon], [180, 90 - _math.epsilon]]).extentMinor([[-180, -80 - _math.epsilon], [180, 80 + _math.epsilon]]);
}

function graticule10() {
  return graticule()();
}
},{"d3-array":"../node_modules/d3-array/src/index.js","./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/interpolate.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _math = require("./math");

function _default(a, b) {
  var x0 = a[0] * _math.radians,
      y0 = a[1] * _math.radians,
      x1 = b[0] * _math.radians,
      y1 = b[1] * _math.radians,
      cy0 = (0, _math.cos)(y0),
      sy0 = (0, _math.sin)(y0),
      cy1 = (0, _math.cos)(y1),
      sy1 = (0, _math.sin)(y1),
      kx0 = cy0 * (0, _math.cos)(x0),
      ky0 = cy0 * (0, _math.sin)(x0),
      kx1 = cy1 * (0, _math.cos)(x1),
      ky1 = cy1 * (0, _math.sin)(x1),
      d = 2 * (0, _math.asin)((0, _math.sqrt)((0, _math.haversin)(y1 - y0) + cy0 * cy1 * (0, _math.haversin)(x1 - x0))),
      k = (0, _math.sin)(d);
  var interpolate = d ? function (t) {
    var B = (0, _math.sin)(t *= d) / k,
        A = (0, _math.sin)(d - t) / k,
        x = A * kx0 + B * kx1,
        y = A * ky0 + B * ky1,
        z = A * sy0 + B * sy1;
    return [(0, _math.atan2)(y, x) * _math.degrees, (0, _math.atan2)(z, (0, _math.sqrt)(x * x + y * y)) * _math.degrees];
  } : function () {
    return [x0 * _math.degrees, y0 * _math.degrees];
  };
  interpolate.distance = d;
  return interpolate;
}
},{"./math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/identity.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return x;
}
},{}],"../node_modules/d3-geo/src/path/area.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _adder = _interopRequireDefault(require("../adder"));

var _math = require("../math");

var _noop = _interopRequireDefault(require("../noop"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var areaSum = (0, _adder.default)(),
    areaRingSum = (0, _adder.default)(),
    x00,
    y00,
    x0,
    y0;
var areaStream = {
  point: _noop.default,
  lineStart: _noop.default,
  lineEnd: _noop.default,
  polygonStart: function () {
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function () {
    areaStream.lineStart = areaStream.lineEnd = areaStream.point = _noop.default;
    areaSum.add((0, _math.abs)(areaRingSum));
    areaRingSum.reset();
  },
  result: function () {
    var area = areaSum / 2;
    areaSum.reset();
    return area;
  }
};

function areaRingStart() {
  areaStream.point = areaPointFirst;
}

function areaPointFirst(x, y) {
  areaStream.point = areaPoint;
  x00 = x0 = x, y00 = y0 = y;
}

function areaPoint(x, y) {
  areaRingSum.add(y0 * x - x0 * y);
  x0 = x, y0 = y;
}

function areaRingEnd() {
  areaPoint(x00, y00);
}

var _default = areaStream;
exports.default = _default;
},{"../adder":"../node_modules/d3-geo/src/adder.js","../math":"../node_modules/d3-geo/src/math.js","../noop":"../node_modules/d3-geo/src/noop.js"}],"../node_modules/d3-geo/src/path/bounds.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _noop = _interopRequireDefault(require("../noop"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var x0 = Infinity,
    y0 = x0,
    x1 = -x0,
    y1 = x1;
var boundsStream = {
  point: boundsPoint,
  lineStart: _noop.default,
  lineEnd: _noop.default,
  polygonStart: _noop.default,
  polygonEnd: _noop.default,
  result: function () {
    var bounds = [[x0, y0], [x1, y1]];
    x1 = y1 = -(y0 = x0 = Infinity);
    return bounds;
  }
};

function boundsPoint(x, y) {
  if (x < x0) x0 = x;
  if (x > x1) x1 = x;
  if (y < y0) y0 = y;
  if (y > y1) y1 = y;
}

var _default = boundsStream;
exports.default = _default;
},{"../noop":"../node_modules/d3-geo/src/noop.js"}],"../node_modules/d3-geo/src/path/centroid.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _math = require("../math");

// TODO Enforce positive area for exterior, negative area for interior?
var X0 = 0,
    Y0 = 0,
    Z0 = 0,
    X1 = 0,
    Y1 = 0,
    Z1 = 0,
    X2 = 0,
    Y2 = 0,
    Z2 = 0,
    x00,
    y00,
    x0,
    y0;
var centroidStream = {
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function () {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function () {
    centroidStream.point = centroidPoint;
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  },
  result: function () {
    var centroid = Z2 ? [X2 / Z2, Y2 / Z2] : Z1 ? [X1 / Z1, Y1 / Z1] : Z0 ? [X0 / Z0, Y0 / Z0] : [NaN, NaN];
    X0 = Y0 = Z0 = X1 = Y1 = Z1 = X2 = Y2 = Z2 = 0;
    return centroid;
  }
};

function centroidPoint(x, y) {
  X0 += x;
  Y0 += y;
  ++Z0;
}

function centroidLineStart() {
  centroidStream.point = centroidPointFirstLine;
}

function centroidPointFirstLine(x, y) {
  centroidStream.point = centroidPointLine;
  centroidPoint(x0 = x, y0 = y);
}

function centroidPointLine(x, y) {
  var dx = x - x0,
      dy = y - y0,
      z = (0, _math.sqrt)(dx * dx + dy * dy);
  X1 += z * (x0 + x) / 2;
  Y1 += z * (y0 + y) / 2;
  Z1 += z;
  centroidPoint(x0 = x, y0 = y);
}

function centroidLineEnd() {
  centroidStream.point = centroidPoint;
}

function centroidRingStart() {
  centroidStream.point = centroidPointFirstRing;
}

function centroidRingEnd() {
  centroidPointRing(x00, y00);
}

function centroidPointFirstRing(x, y) {
  centroidStream.point = centroidPointRing;
  centroidPoint(x00 = x0 = x, y00 = y0 = y);
}

function centroidPointRing(x, y) {
  var dx = x - x0,
      dy = y - y0,
      z = (0, _math.sqrt)(dx * dx + dy * dy);
  X1 += z * (x0 + x) / 2;
  Y1 += z * (y0 + y) / 2;
  Z1 += z;
  z = y0 * x - x0 * y;
  X2 += z * (x0 + x);
  Y2 += z * (y0 + y);
  Z2 += z * 3;
  centroidPoint(x0 = x, y0 = y);
}

var _default = centroidStream;
exports.default = _default;
},{"../math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/path/context.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = PathContext;

var _math = require("../math");

var _noop = _interopRequireDefault(require("../noop"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function PathContext(context) {
  this._context = context;
}

PathContext.prototype = {
  _radius: 4.5,
  pointRadius: function (_) {
    return this._radius = _, this;
  },
  polygonStart: function () {
    this._line = 0;
  },
  polygonEnd: function () {
    this._line = NaN;
  },
  lineStart: function () {
    this._point = 0;
  },
  lineEnd: function () {
    if (this._line === 0) this._context.closePath();
    this._point = NaN;
  },
  point: function (x, y) {
    switch (this._point) {
      case 0:
        {
          this._context.moveTo(x, y);

          this._point = 1;
          break;
        }

      case 1:
        {
          this._context.lineTo(x, y);

          break;
        }

      default:
        {
          this._context.moveTo(x + this._radius, y);

          this._context.arc(x, y, this._radius, 0, _math.tau);

          break;
        }
    }
  },
  result: _noop.default
};
},{"../math":"../node_modules/d3-geo/src/math.js","../noop":"../node_modules/d3-geo/src/noop.js"}],"../node_modules/d3-geo/src/path/measure.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _adder = _interopRequireDefault(require("../adder"));

var _math = require("../math");

var _noop = _interopRequireDefault(require("../noop"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var lengthSum = (0, _adder.default)(),
    lengthRing,
    x00,
    y00,
    x0,
    y0;
var lengthStream = {
  point: _noop.default,
  lineStart: function () {
    lengthStream.point = lengthPointFirst;
  },
  lineEnd: function () {
    if (lengthRing) lengthPoint(x00, y00);
    lengthStream.point = _noop.default;
  },
  polygonStart: function () {
    lengthRing = true;
  },
  polygonEnd: function () {
    lengthRing = null;
  },
  result: function () {
    var length = +lengthSum;
    lengthSum.reset();
    return length;
  }
};

function lengthPointFirst(x, y) {
  lengthStream.point = lengthPoint;
  x00 = x0 = x, y00 = y0 = y;
}

function lengthPoint(x, y) {
  x0 -= x, y0 -= y;
  lengthSum.add((0, _math.sqrt)(x0 * x0 + y0 * y0));
  x0 = x, y0 = y;
}

var _default = lengthStream;
exports.default = _default;
},{"../adder":"../node_modules/d3-geo/src/adder.js","../math":"../node_modules/d3-geo/src/math.js","../noop":"../node_modules/d3-geo/src/noop.js"}],"../node_modules/d3-geo/src/path/string.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = PathString;

function PathString() {
  this._string = [];
}

PathString.prototype = {
  _radius: 4.5,
  _circle: circle(4.5),
  pointRadius: function (_) {
    if ((_ = +_) !== this._radius) this._radius = _, this._circle = null;
    return this;
  },
  polygonStart: function () {
    this._line = 0;
  },
  polygonEnd: function () {
    this._line = NaN;
  },
  lineStart: function () {
    this._point = 0;
  },
  lineEnd: function () {
    if (this._line === 0) this._string.push("Z");
    this._point = NaN;
  },
  point: function (x, y) {
    switch (this._point) {
      case 0:
        {
          this._string.push("M", x, ",", y);

          this._point = 1;
          break;
        }

      case 1:
        {
          this._string.push("L", x, ",", y);

          break;
        }

      default:
        {
          if (this._circle == null) this._circle = circle(this._radius);

          this._string.push("M", x, ",", y, this._circle);

          break;
        }
    }
  },
  result: function () {
    if (this._string.length) {
      var result = this._string.join("");

      this._string = [];
      return result;
    } else {
      return null;
    }
  }
};

function circle(radius) {
  return "m0," + radius + "a" + radius + "," + radius + " 0 1,1 0," + -2 * radius + "a" + radius + "," + radius + " 0 1,1 0," + 2 * radius + "z";
}
},{}],"../node_modules/d3-geo/src/path/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _identity = _interopRequireDefault(require("../identity"));

var _stream = _interopRequireDefault(require("../stream"));

var _area = _interopRequireDefault(require("./area"));

var _bounds = _interopRequireDefault(require("./bounds"));

var _centroid = _interopRequireDefault(require("./centroid"));

var _context = _interopRequireDefault(require("./context"));

var _measure = _interopRequireDefault(require("./measure"));

var _string = _interopRequireDefault(require("./string"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(projection, context) {
  var pointRadius = 4.5,
      projectionStream,
      contextStream;

  function path(object) {
    if (object) {
      if (typeof pointRadius === "function") contextStream.pointRadius(+pointRadius.apply(this, arguments));
      (0, _stream.default)(object, projectionStream(contextStream));
    }

    return contextStream.result();
  }

  path.area = function (object) {
    (0, _stream.default)(object, projectionStream(_area.default));
    return _area.default.result();
  };

  path.measure = function (object) {
    (0, _stream.default)(object, projectionStream(_measure.default));
    return _measure.default.result();
  };

  path.bounds = function (object) {
    (0, _stream.default)(object, projectionStream(_bounds.default));
    return _bounds.default.result();
  };

  path.centroid = function (object) {
    (0, _stream.default)(object, projectionStream(_centroid.default));
    return _centroid.default.result();
  };

  path.projection = function (_) {
    return arguments.length ? (projectionStream = _ == null ? (projection = null, _identity.default) : (projection = _).stream, path) : projection;
  };

  path.context = function (_) {
    if (!arguments.length) return context;
    contextStream = _ == null ? (context = null, new _string.default()) : new _context.default(context = _);
    if (typeof pointRadius !== "function") contextStream.pointRadius(pointRadius);
    return path;
  };

  path.pointRadius = function (_) {
    if (!arguments.length) return pointRadius;
    pointRadius = typeof _ === "function" ? _ : (contextStream.pointRadius(+_), +_);
    return path;
  };

  return path.projection(projection).context(context);
}
},{"../identity":"../node_modules/d3-geo/src/identity.js","../stream":"../node_modules/d3-geo/src/stream.js","./area":"../node_modules/d3-geo/src/path/area.js","./bounds":"../node_modules/d3-geo/src/path/bounds.js","./centroid":"../node_modules/d3-geo/src/path/centroid.js","./context":"../node_modules/d3-geo/src/path/context.js","./measure":"../node_modules/d3-geo/src/path/measure.js","./string":"../node_modules/d3-geo/src/path/string.js"}],"../node_modules/d3-geo/src/transform.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.transformer = transformer;

function _default(methods) {
  return {
    stream: transformer(methods)
  };
}

function transformer(methods) {
  return function (stream) {
    var s = new TransformStream();

    for (var key in methods) s[key] = methods[key];

    s.stream = stream;
    return s;
  };
}

function TransformStream() {}

TransformStream.prototype = {
  constructor: TransformStream,
  point: function (x, y) {
    this.stream.point(x, y);
  },
  sphere: function () {
    this.stream.sphere();
  },
  lineStart: function () {
    this.stream.lineStart();
  },
  lineEnd: function () {
    this.stream.lineEnd();
  },
  polygonStart: function () {
    this.stream.polygonStart();
  },
  polygonEnd: function () {
    this.stream.polygonEnd();
  }
};
},{}],"../node_modules/d3-geo/src/projection/fit.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.fitExtent = fitExtent;
exports.fitSize = fitSize;
exports.fitWidth = fitWidth;
exports.fitHeight = fitHeight;

var _stream = _interopRequireDefault(require("../stream"));

var _bounds = _interopRequireDefault(require("../path/bounds"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function fit(projection, fitBounds, object) {
  var clip = projection.clipExtent && projection.clipExtent();
  projection.scale(150).translate([0, 0]);
  if (clip != null) projection.clipExtent(null);
  (0, _stream.default)(object, projection.stream(_bounds.default));
  fitBounds(_bounds.default.result());
  if (clip != null) projection.clipExtent(clip);
  return projection;
}

function fitExtent(projection, extent, object) {
  return fit(projection, function (b) {
    var w = extent[1][0] - extent[0][0],
        h = extent[1][1] - extent[0][1],
        k = Math.min(w / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])),
        x = +extent[0][0] + (w - k * (b[1][0] + b[0][0])) / 2,
        y = +extent[0][1] + (h - k * (b[1][1] + b[0][1])) / 2;
    projection.scale(150 * k).translate([x, y]);
  }, object);
}

function fitSize(projection, size, object) {
  return fitExtent(projection, [[0, 0], size], object);
}

function fitWidth(projection, width, object) {
  return fit(projection, function (b) {
    var w = +width,
        k = w / (b[1][0] - b[0][0]),
        x = (w - k * (b[1][0] + b[0][0])) / 2,
        y = -k * b[0][1];
    projection.scale(150 * k).translate([x, y]);
  }, object);
}

function fitHeight(projection, height, object) {
  return fit(projection, function (b) {
    var h = +height,
        k = h / (b[1][1] - b[0][1]),
        x = -k * b[0][0],
        y = (h - k * (b[1][1] + b[0][1])) / 2;
    projection.scale(150 * k).translate([x, y]);
  }, object);
}
},{"../stream":"../node_modules/d3-geo/src/stream.js","../path/bounds":"../node_modules/d3-geo/src/path/bounds.js"}],"../node_modules/d3-geo/src/projection/resample.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _cartesian = require("../cartesian");

var _math = require("../math");

var _transform = require("../transform");

var maxDepth = 16,
    // maximum depth of subdivision
cosMinDistance = (0, _math.cos)(30 * _math.radians); // cos(minimum angular distance)

function _default(project, delta2) {
  return +delta2 ? resample(project, delta2) : resampleNone(project);
}

function resampleNone(project) {
  return (0, _transform.transformer)({
    point: function (x, y) {
      x = project(x, y);
      this.stream.point(x[0], x[1]);
    }
  });
}

function resample(project, delta2) {
  function resampleLineTo(x0, y0, lambda0, a0, b0, c0, x1, y1, lambda1, a1, b1, c1, depth, stream) {
    var dx = x1 - x0,
        dy = y1 - y0,
        d2 = dx * dx + dy * dy;

    if (d2 > 4 * delta2 && depth--) {
      var a = a0 + a1,
          b = b0 + b1,
          c = c0 + c1,
          m = (0, _math.sqrt)(a * a + b * b + c * c),
          phi2 = (0, _math.asin)(c /= m),
          lambda2 = (0, _math.abs)((0, _math.abs)(c) - 1) < _math.epsilon || (0, _math.abs)(lambda0 - lambda1) < _math.epsilon ? (lambda0 + lambda1) / 2 : (0, _math.atan2)(b, a),
          p = project(lambda2, phi2),
          x2 = p[0],
          y2 = p[1],
          dx2 = x2 - x0,
          dy2 = y2 - y0,
          dz = dy * dx2 - dx * dy2;

      if (dz * dz / d2 > delta2 // perpendicular projected distance
      || (0, _math.abs)((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 // midpoint close to an end
      || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) {
        // angular distance
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x2, y2, lambda2, a /= m, b /= m, c, depth, stream);
        stream.point(x2, y2);
        resampleLineTo(x2, y2, lambda2, a, b, c, x1, y1, lambda1, a1, b1, c1, depth, stream);
      }
    }
  }

  return function (stream) {
    var lambda00, x00, y00, a00, b00, c00, // first point
    lambda0, x0, y0, a0, b0, c0; // previous point

    var resampleStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function () {
        stream.polygonStart();
        resampleStream.lineStart = ringStart;
      },
      polygonEnd: function () {
        stream.polygonEnd();
        resampleStream.lineStart = lineStart;
      }
    };

    function point(x, y) {
      x = project(x, y);
      stream.point(x[0], x[1]);
    }

    function lineStart() {
      x0 = NaN;
      resampleStream.point = linePoint;
      stream.lineStart();
    }

    function linePoint(lambda, phi) {
      var c = (0, _cartesian.cartesian)([lambda, phi]),
          p = project(lambda, phi);
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x0 = p[0], y0 = p[1], lambda0 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
      stream.point(x0, y0);
    }

    function lineEnd() {
      resampleStream.point = point;
      stream.lineEnd();
    }

    function ringStart() {
      lineStart();
      resampleStream.point = ringPoint;
      resampleStream.lineEnd = ringEnd;
    }

    function ringPoint(lambda, phi) {
      linePoint(lambda00 = lambda, phi), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
      resampleStream.point = linePoint;
    }

    function ringEnd() {
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x00, y00, lambda00, a00, b00, c00, maxDepth, stream);
      resampleStream.lineEnd = lineEnd;
      lineEnd();
    }

    return resampleStream;
  };
}
},{"../cartesian":"../node_modules/d3-geo/src/cartesian.js","../math":"../node_modules/d3-geo/src/math.js","../transform":"../node_modules/d3-geo/src/transform.js"}],"../node_modules/d3-geo/src/projection/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = projection;
exports.projectionMutator = projectionMutator;

var _antimeridian = _interopRequireDefault(require("../clip/antimeridian"));

var _circle = _interopRequireDefault(require("../clip/circle"));

var _rectangle = _interopRequireDefault(require("../clip/rectangle"));

var _compose = _interopRequireDefault(require("../compose"));

var _identity = _interopRequireDefault(require("../identity"));

var _math = require("../math");

var _rotation = require("../rotation");

var _transform = require("../transform");

var _fit = require("./fit");

var _resample = _interopRequireDefault(require("./resample"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var transformRadians = (0, _transform.transformer)({
  point: function (x, y) {
    this.stream.point(x * _math.radians, y * _math.radians);
  }
});

function transformRotate(rotate) {
  return (0, _transform.transformer)({
    point: function (x, y) {
      var r = rotate(x, y);
      return this.stream.point(r[0], r[1]);
    }
  });
}

function scaleTranslate(k, dx, dy) {
  function transform(x, y) {
    return [dx + k * x, dy - k * y];
  }

  transform.invert = function (x, y) {
    return [(x - dx) / k, (dy - y) / k];
  };

  return transform;
}

function scaleTranslateRotate(k, dx, dy, alpha) {
  var cosAlpha = (0, _math.cos)(alpha),
      sinAlpha = (0, _math.sin)(alpha),
      a = cosAlpha * k,
      b = sinAlpha * k,
      ai = cosAlpha / k,
      bi = sinAlpha / k,
      ci = (sinAlpha * dy - cosAlpha * dx) / k,
      fi = (sinAlpha * dx + cosAlpha * dy) / k;

  function transform(x, y) {
    return [a * x - b * y + dx, dy - b * x - a * y];
  }

  transform.invert = function (x, y) {
    return [ai * x - bi * y + ci, fi - bi * x - ai * y];
  };

  return transform;
}

function projection(project) {
  return projectionMutator(function () {
    return project;
  })();
}

function projectionMutator(projectAt) {
  var project,
      k = 150,
      // scale
  x = 480,
      y = 250,
      // translate
  lambda = 0,
      phi = 0,
      // center
  deltaLambda = 0,
      deltaPhi = 0,
      deltaGamma = 0,
      rotate,
      // pre-rotate
  alpha = 0,
      // post-rotate
  theta = null,
      preclip = _antimeridian.default,
      // pre-clip angle
  x0 = null,
      y0,
      x1,
      y1,
      postclip = _identity.default,
      // post-clip extent
  delta2 = 0.5,
      // precision
  projectResample,
      projectTransform,
      projectRotateTransform,
      cache,
      cacheStream;

  function projection(point) {
    return projectRotateTransform(point[0] * _math.radians, point[1] * _math.radians);
  }

  function invert(point) {
    point = projectRotateTransform.invert(point[0], point[1]);
    return point && [point[0] * _math.degrees, point[1] * _math.degrees];
  }

  projection.stream = function (stream) {
    return cache && cacheStream === stream ? cache : cache = transformRadians(transformRotate(rotate)(preclip(projectResample(postclip(cacheStream = stream)))));
  };

  projection.preclip = function (_) {
    return arguments.length ? (preclip = _, theta = undefined, reset()) : preclip;
  };

  projection.postclip = function (_) {
    return arguments.length ? (postclip = _, x0 = y0 = x1 = y1 = null, reset()) : postclip;
  };

  projection.clipAngle = function (_) {
    return arguments.length ? (preclip = +_ ? (0, _circle.default)(theta = _ * _math.radians) : (theta = null, _antimeridian.default), reset()) : theta * _math.degrees;
  };

  projection.clipExtent = function (_) {
    return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, _identity.default) : (0, _rectangle.default)(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  projection.scale = function (_) {
    return arguments.length ? (k = +_, recenter()) : k;
  };

  projection.translate = function (_) {
    return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
  };

  projection.center = function (_) {
    return arguments.length ? (lambda = _[0] % 360 * _math.radians, phi = _[1] % 360 * _math.radians, recenter()) : [lambda * _math.degrees, phi * _math.degrees];
  };

  projection.rotate = function (_) {
    return arguments.length ? (deltaLambda = _[0] % 360 * _math.radians, deltaPhi = _[1] % 360 * _math.radians, deltaGamma = _.length > 2 ? _[2] % 360 * _math.radians : 0, recenter()) : [deltaLambda * _math.degrees, deltaPhi * _math.degrees, deltaGamma * _math.degrees];
  };

  projection.angle = function (_) {
    return arguments.length ? (alpha = _ % 360 * _math.radians, recenter()) : alpha * _math.degrees;
  };

  projection.precision = function (_) {
    return arguments.length ? (projectResample = (0, _resample.default)(projectTransform, delta2 = _ * _), reset()) : (0, _math.sqrt)(delta2);
  };

  projection.fitExtent = function (extent, object) {
    return (0, _fit.fitExtent)(projection, extent, object);
  };

  projection.fitSize = function (size, object) {
    return (0, _fit.fitSize)(projection, size, object);
  };

  projection.fitWidth = function (width, object) {
    return (0, _fit.fitWidth)(projection, width, object);
  };

  projection.fitHeight = function (height, object) {
    return (0, _fit.fitHeight)(projection, height, object);
  };

  function recenter() {
    var center = scaleTranslateRotate(k, 0, 0, alpha).apply(null, project(lambda, phi)),
        transform = (alpha ? scaleTranslateRotate : scaleTranslate)(k, x - center[0], y - center[1], alpha);
    rotate = (0, _rotation.rotateRadians)(deltaLambda, deltaPhi, deltaGamma);
    projectTransform = (0, _compose.default)(project, transform);
    projectRotateTransform = (0, _compose.default)(rotate, projectTransform);
    projectResample = (0, _resample.default)(projectTransform, delta2);
    return reset();
  }

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return function () {
    project = projectAt.apply(this, arguments);
    projection.invert = project.invert && invert;
    return recenter();
  };
}
},{"../clip/antimeridian":"../node_modules/d3-geo/src/clip/antimeridian.js","../clip/circle":"../node_modules/d3-geo/src/clip/circle.js","../clip/rectangle":"../node_modules/d3-geo/src/clip/rectangle.js","../compose":"../node_modules/d3-geo/src/compose.js","../identity":"../node_modules/d3-geo/src/identity.js","../math":"../node_modules/d3-geo/src/math.js","../rotation":"../node_modules/d3-geo/src/rotation.js","../transform":"../node_modules/d3-geo/src/transform.js","./fit":"../node_modules/d3-geo/src/projection/fit.js","./resample":"../node_modules/d3-geo/src/projection/resample.js"}],"../node_modules/d3-geo/src/projection/conic.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.conicProjection = conicProjection;

var _math = require("../math");

var _index = require("./index");

function conicProjection(projectAt) {
  var phi0 = 0,
      phi1 = _math.pi / 3,
      m = (0, _index.projectionMutator)(projectAt),
      p = m(phi0, phi1);

  p.parallels = function (_) {
    return arguments.length ? m(phi0 = _[0] * _math.radians, phi1 = _[1] * _math.radians) : [phi0 * _math.degrees, phi1 * _math.degrees];
  };

  return p;
}
},{"../math":"../node_modules/d3-geo/src/math.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/cylindricalEqualArea.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.cylindricalEqualAreaRaw = cylindricalEqualAreaRaw;

var _math = require("../math");

function cylindricalEqualAreaRaw(phi0) {
  var cosPhi0 = (0, _math.cos)(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, (0, _math.sin)(phi) / cosPhi0];
  }

  forward.invert = function (x, y) {
    return [x / cosPhi0, (0, _math.asin)(y * cosPhi0)];
  };

  return forward;
}
},{"../math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/projection/conicEqualArea.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.conicEqualAreaRaw = conicEqualAreaRaw;
exports.default = _default;

var _math = require("../math");

var _conic = require("./conic");

var _cylindricalEqualArea = require("./cylindricalEqualArea");

function conicEqualAreaRaw(y0, y1) {
  var sy0 = (0, _math.sin)(y0),
      n = (sy0 + (0, _math.sin)(y1)) / 2; // Are the parallels symmetrical around the Equator?

  if ((0, _math.abs)(n) < _math.epsilon) return (0, _cylindricalEqualArea.cylindricalEqualAreaRaw)(y0);
  var c = 1 + sy0 * (2 * n - sy0),
      r0 = (0, _math.sqrt)(c) / n;

  function project(x, y) {
    var r = (0, _math.sqrt)(c - 2 * n * (0, _math.sin)(y)) / n;
    return [r * (0, _math.sin)(x *= n), r0 - r * (0, _math.cos)(x)];
  }

  project.invert = function (x, y) {
    var r0y = r0 - y;
    return [(0, _math.atan2)(x, (0, _math.abs)(r0y)) / n * (0, _math.sign)(r0y), (0, _math.asin)((c - (x * x + r0y * r0y) * n * n) / (2 * n))];
  };

  return project;
}

function _default() {
  return (0, _conic.conicProjection)(conicEqualAreaRaw).scale(155.424).center([0, 33.6442]);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./conic":"../node_modules/d3-geo/src/projection/conic.js","./cylindricalEqualArea":"../node_modules/d3-geo/src/projection/cylindricalEqualArea.js"}],"../node_modules/d3-geo/src/projection/albers.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _conicEqualArea = _interopRequireDefault(require("./conicEqualArea"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  return (0, _conicEqualArea.default)().parallels([29.5, 45.5]).scale(1070).translate([480, 250]).rotate([96, 0]).center([-0.6, 38.7]);
}
},{"./conicEqualArea":"../node_modules/d3-geo/src/projection/conicEqualArea.js"}],"../node_modules/d3-geo/src/projection/albersUsa.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _math = require("../math");

var _albers = _interopRequireDefault(require("./albers"));

var _conicEqualArea = _interopRequireDefault(require("./conicEqualArea"));

var _fit = require("./fit");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

// The projections must have mutually exclusive clip regions on the sphere,
// as this will avoid emitting interleaving lines and polygons.
function multiplex(streams) {
  var n = streams.length;
  return {
    point: function (x, y) {
      var i = -1;

      while (++i < n) streams[i].point(x, y);
    },
    sphere: function () {
      var i = -1;

      while (++i < n) streams[i].sphere();
    },
    lineStart: function () {
      var i = -1;

      while (++i < n) streams[i].lineStart();
    },
    lineEnd: function () {
      var i = -1;

      while (++i < n) streams[i].lineEnd();
    },
    polygonStart: function () {
      var i = -1;

      while (++i < n) streams[i].polygonStart();
    },
    polygonEnd: function () {
      var i = -1;

      while (++i < n) streams[i].polygonEnd();
    }
  };
} // A composite projection for the United States, configured by default for
// 960×500. The projection also works quite well at 960×600 if you change the
// scale to 1285 and adjust the translate accordingly. The set of standard
// parallels for each region comes from USGS, which is published here:
// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers


function _default() {
  var cache,
      cacheStream,
      lower48 = (0, _albers.default)(),
      lower48Point,
      alaska = (0, _conicEqualArea.default)().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]),
      alaskaPoint,
      // EPSG:3338
  hawaii = (0, _conicEqualArea.default)().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]),
      hawaiiPoint,
      // ESRI:102007
  point,
      pointStream = {
    point: function (x, y) {
      point = [x, y];
    }
  };

  function albersUsa(coordinates) {
    var x = coordinates[0],
        y = coordinates[1];
    return point = null, (lower48Point.point(x, y), point) || (alaskaPoint.point(x, y), point) || (hawaiiPoint.point(x, y), point);
  }

  albersUsa.invert = function (coordinates) {
    var k = lower48.scale(),
        t = lower48.translate(),
        x = (coordinates[0] - t[0]) / k,
        y = (coordinates[1] - t[1]) / k;
    return (y >= 0.120 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii : lower48).invert(coordinates);
  };

  albersUsa.stream = function (stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex([lower48.stream(cacheStream = stream), alaska.stream(stream), hawaii.stream(stream)]);
  };

  albersUsa.precision = function (_) {
    if (!arguments.length) return lower48.precision();
    lower48.precision(_), alaska.precision(_), hawaii.precision(_);
    return reset();
  };

  albersUsa.scale = function (_) {
    if (!arguments.length) return lower48.scale();
    lower48.scale(_), alaska.scale(_ * 0.35), hawaii.scale(_);
    return albersUsa.translate(lower48.translate());
  };

  albersUsa.translate = function (_) {
    if (!arguments.length) return lower48.translate();
    var k = lower48.scale(),
        x = +_[0],
        y = +_[1];
    lower48Point = lower48.translate(_).clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]]).stream(pointStream);
    alaskaPoint = alaska.translate([x - 0.307 * k, y + 0.201 * k]).clipExtent([[x - 0.425 * k + _math.epsilon, y + 0.120 * k + _math.epsilon], [x - 0.214 * k - _math.epsilon, y + 0.234 * k - _math.epsilon]]).stream(pointStream);
    hawaiiPoint = hawaii.translate([x - 0.205 * k, y + 0.212 * k]).clipExtent([[x - 0.214 * k + _math.epsilon, y + 0.166 * k + _math.epsilon], [x - 0.115 * k - _math.epsilon, y + 0.234 * k - _math.epsilon]]).stream(pointStream);
    return reset();
  };

  albersUsa.fitExtent = function (extent, object) {
    return (0, _fit.fitExtent)(albersUsa, extent, object);
  };

  albersUsa.fitSize = function (size, object) {
    return (0, _fit.fitSize)(albersUsa, size, object);
  };

  albersUsa.fitWidth = function (width, object) {
    return (0, _fit.fitWidth)(albersUsa, width, object);
  };

  albersUsa.fitHeight = function (height, object) {
    return (0, _fit.fitHeight)(albersUsa, height, object);
  };

  function reset() {
    cache = cacheStream = null;
    return albersUsa;
  }

  return albersUsa.scale(1070);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./albers":"../node_modules/d3-geo/src/projection/albers.js","./conicEqualArea":"../node_modules/d3-geo/src/projection/conicEqualArea.js","./fit":"../node_modules/d3-geo/src/projection/fit.js"}],"../node_modules/d3-geo/src/projection/azimuthal.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.azimuthalRaw = azimuthalRaw;
exports.azimuthalInvert = azimuthalInvert;

var _math = require("../math");

function azimuthalRaw(scale) {
  return function (x, y) {
    var cx = (0, _math.cos)(x),
        cy = (0, _math.cos)(y),
        k = scale(cx * cy);
    return [k * cy * (0, _math.sin)(x), k * (0, _math.sin)(y)];
  };
}

function azimuthalInvert(angle) {
  return function (x, y) {
    var z = (0, _math.sqrt)(x * x + y * y),
        c = angle(z),
        sc = (0, _math.sin)(c),
        cc = (0, _math.cos)(c);
    return [(0, _math.atan2)(x * sc, z * cc), (0, _math.asin)(z && y * sc / z)];
  };
}
},{"../math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/projection/azimuthalEqualArea.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.azimuthalEqualAreaRaw = void 0;

var _math = require("../math");

var _azimuthal = require("./azimuthal");

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var azimuthalEqualAreaRaw = (0, _azimuthal.azimuthalRaw)(function (cxcy) {
  return (0, _math.sqrt)(2 / (1 + cxcy));
});
exports.azimuthalEqualAreaRaw = azimuthalEqualAreaRaw;
azimuthalEqualAreaRaw.invert = (0, _azimuthal.azimuthalInvert)(function (z) {
  return 2 * (0, _math.asin)(z / 2);
});

function _default() {
  return (0, _index.default)(azimuthalEqualAreaRaw).scale(124.75).clipAngle(180 - 1e-3);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./azimuthal":"../node_modules/d3-geo/src/projection/azimuthal.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/azimuthalEquidistant.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.azimuthalEquidistantRaw = void 0;

var _math = require("../math");

var _azimuthal = require("./azimuthal");

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var azimuthalEquidistantRaw = (0, _azimuthal.azimuthalRaw)(function (c) {
  return (c = (0, _math.acos)(c)) && c / (0, _math.sin)(c);
});
exports.azimuthalEquidistantRaw = azimuthalEquidistantRaw;
azimuthalEquidistantRaw.invert = (0, _azimuthal.azimuthalInvert)(function (z) {
  return z;
});

function _default() {
  return (0, _index.default)(azimuthalEquidistantRaw).scale(79.4188).clipAngle(180 - 1e-3);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./azimuthal":"../node_modules/d3-geo/src/projection/azimuthal.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/mercator.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.mercatorRaw = mercatorRaw;
exports.default = _default;
exports.mercatorProjection = mercatorProjection;

var _math = require("../math");

var _rotation = _interopRequireDefault(require("../rotation"));

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function mercatorRaw(lambda, phi) {
  return [lambda, (0, _math.log)((0, _math.tan)((_math.halfPi + phi) / 2))];
}

mercatorRaw.invert = function (x, y) {
  return [x, 2 * (0, _math.atan)((0, _math.exp)(y)) - _math.halfPi];
};

function _default() {
  return mercatorProjection(mercatorRaw).scale(961 / _math.tau);
}

function mercatorProjection(project) {
  var m = (0, _index.default)(project),
      center = m.center,
      scale = m.scale,
      translate = m.translate,
      clipExtent = m.clipExtent,
      x0 = null,
      y0,
      x1,
      y1; // clip extent

  m.scale = function (_) {
    return arguments.length ? (scale(_), reclip()) : scale();
  };

  m.translate = function (_) {
    return arguments.length ? (translate(_), reclip()) : translate();
  };

  m.center = function (_) {
    return arguments.length ? (center(_), reclip()) : center();
  };

  m.clipExtent = function (_) {
    return arguments.length ? (_ == null ? x0 = y0 = x1 = y1 = null : (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reclip()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  function reclip() {
    var k = _math.pi * scale(),
        t = m((0, _rotation.default)(m.rotate()).invert([0, 0]));
    return clipExtent(x0 == null ? [[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]] : project === mercatorRaw ? [[Math.max(t[0] - k, x0), y0], [Math.min(t[0] + k, x1), y1]] : [[x0, Math.max(t[1] - k, y0)], [x1, Math.min(t[1] + k, y1)]]);
  }

  return reclip();
}
},{"../math":"../node_modules/d3-geo/src/math.js","../rotation":"../node_modules/d3-geo/src/rotation.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/conicConformal.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.conicConformalRaw = conicConformalRaw;
exports.default = _default;

var _math = require("../math");

var _conic = require("./conic");

var _mercator = require("./mercator");

function tany(y) {
  return (0, _math.tan)((_math.halfPi + y) / 2);
}

function conicConformalRaw(y0, y1) {
  var cy0 = (0, _math.cos)(y0),
      n = y0 === y1 ? (0, _math.sin)(y0) : (0, _math.log)(cy0 / (0, _math.cos)(y1)) / (0, _math.log)(tany(y1) / tany(y0)),
      f = cy0 * (0, _math.pow)(tany(y0), n) / n;
  if (!n) return _mercator.mercatorRaw;

  function project(x, y) {
    if (f > 0) {
      if (y < -_math.halfPi + _math.epsilon) y = -_math.halfPi + _math.epsilon;
    } else {
      if (y > _math.halfPi - _math.epsilon) y = _math.halfPi - _math.epsilon;
    }

    var r = f / (0, _math.pow)(tany(y), n);
    return [r * (0, _math.sin)(n * x), f - r * (0, _math.cos)(n * x)];
  }

  project.invert = function (x, y) {
    var fy = f - y,
        r = (0, _math.sign)(n) * (0, _math.sqrt)(x * x + fy * fy);
    return [(0, _math.atan2)(x, (0, _math.abs)(fy)) / n * (0, _math.sign)(fy), 2 * (0, _math.atan)((0, _math.pow)(f / r, 1 / n)) - _math.halfPi];
  };

  return project;
}

function _default() {
  return (0, _conic.conicProjection)(conicConformalRaw).scale(109.5).parallels([30, 30]);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./conic":"../node_modules/d3-geo/src/projection/conic.js","./mercator":"../node_modules/d3-geo/src/projection/mercator.js"}],"../node_modules/d3-geo/src/projection/equirectangular.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.equirectangularRaw = equirectangularRaw;
exports.default = _default;

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function equirectangularRaw(lambda, phi) {
  return [lambda, phi];
}

equirectangularRaw.invert = equirectangularRaw;

function _default() {
  return (0, _index.default)(equirectangularRaw).scale(152.63);
}
},{"./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/conicEquidistant.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.conicEquidistantRaw = conicEquidistantRaw;
exports.default = _default;

var _math = require("../math");

var _conic = require("./conic");

var _equirectangular = require("./equirectangular");

function conicEquidistantRaw(y0, y1) {
  var cy0 = (0, _math.cos)(y0),
      n = y0 === y1 ? (0, _math.sin)(y0) : (cy0 - (0, _math.cos)(y1)) / (y1 - y0),
      g = cy0 / n + y0;
  if ((0, _math.abs)(n) < _math.epsilon) return _equirectangular.equirectangularRaw;

  function project(x, y) {
    var gy = g - y,
        nx = n * x;
    return [gy * (0, _math.sin)(nx), g - gy * (0, _math.cos)(nx)];
  }

  project.invert = function (x, y) {
    var gy = g - y;
    return [(0, _math.atan2)(x, (0, _math.abs)(gy)) / n * (0, _math.sign)(gy), g - (0, _math.sign)(n) * (0, _math.sqrt)(x * x + gy * gy)];
  };

  return project;
}

function _default() {
  return (0, _conic.conicProjection)(conicEquidistantRaw).scale(131.154).center([0, 13.9389]);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./conic":"../node_modules/d3-geo/src/projection/conic.js","./equirectangular":"../node_modules/d3-geo/src/projection/equirectangular.js"}],"../node_modules/d3-geo/src/projection/equalEarth.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.equalEarthRaw = equalEarthRaw;
exports.default = _default;

var _index = _interopRequireDefault(require("./index.js"));

var _math = require("../math.js");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var A1 = 1.340264,
    A2 = -0.081106,
    A3 = 0.000893,
    A4 = 0.003796,
    M = (0, _math.sqrt)(3) / 2,
    iterations = 12;

function equalEarthRaw(lambda, phi) {
  var l = (0, _math.asin)(M * (0, _math.sin)(phi)),
      l2 = l * l,
      l6 = l2 * l2 * l2;
  return [lambda * (0, _math.cos)(l) / (M * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2))), l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2))];
}

equalEarthRaw.invert = function (x, y) {
  var l = y,
      l2 = l * l,
      l6 = l2 * l2 * l2;

  for (var i = 0, delta, fy, fpy; i < iterations; ++i) {
    fy = l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2)) - y;
    fpy = A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2);
    l -= delta = fy / fpy, l2 = l * l, l6 = l2 * l2 * l2;
    if ((0, _math.abs)(delta) < _math.epsilon2) break;
  }

  return [M * x * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2)) / (0, _math.cos)(l), (0, _math.asin)((0, _math.sin)(l) / M)];
};

function _default() {
  return (0, _index.default)(equalEarthRaw).scale(177.158);
}
},{"./index.js":"../node_modules/d3-geo/src/projection/index.js","../math.js":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/projection/gnomonic.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.gnomonicRaw = gnomonicRaw;
exports.default = _default;

var _math = require("../math");

var _azimuthal = require("./azimuthal");

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function gnomonicRaw(x, y) {
  var cy = (0, _math.cos)(y),
      k = (0, _math.cos)(x) * cy;
  return [cy * (0, _math.sin)(x) / k, (0, _math.sin)(y) / k];
}

gnomonicRaw.invert = (0, _azimuthal.azimuthalInvert)(_math.atan);

function _default() {
  return (0, _index.default)(gnomonicRaw).scale(144.049).clipAngle(60);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./azimuthal":"../node_modules/d3-geo/src/projection/azimuthal.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/identity.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _rectangle = _interopRequireDefault(require("../clip/rectangle"));

var _identity = _interopRequireDefault(require("../identity"));

var _transform = require("../transform");

var _fit = require("./fit");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function scaleTranslate(kx, ky, tx, ty) {
  return kx === 1 && ky === 1 && tx === 0 && ty === 0 ? _identity.default : (0, _transform.transformer)({
    point: function (x, y) {
      this.stream.point(x * kx + tx, y * ky + ty);
    }
  });
}

function _default() {
  var k = 1,
      tx = 0,
      ty = 0,
      sx = 1,
      sy = 1,
      transform = _identity.default,
      // scale, translate and reflect
  x0 = null,
      y0,
      x1,
      y1,
      // clip extent
  postclip = _identity.default,
      cache,
      cacheStream,
      projection;

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return projection = {
    stream: function (stream) {
      return cache && cacheStream === stream ? cache : cache = transform(postclip(cacheStream = stream));
    },
    postclip: function (_) {
      return arguments.length ? (postclip = _, x0 = y0 = x1 = y1 = null, reset()) : postclip;
    },
    clipExtent: function (_) {
      return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, _identity.default) : (0, _rectangle.default)(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
    },
    scale: function (_) {
      return arguments.length ? (transform = scaleTranslate((k = +_) * sx, k * sy, tx, ty), reset()) : k;
    },
    translate: function (_) {
      return arguments.length ? (transform = scaleTranslate(k * sx, k * sy, tx = +_[0], ty = +_[1]), reset()) : [tx, ty];
    },
    reflectX: function (_) {
      return arguments.length ? (transform = scaleTranslate(k * (sx = _ ? -1 : 1), k * sy, tx, ty), reset()) : sx < 0;
    },
    reflectY: function (_) {
      return arguments.length ? (transform = scaleTranslate(k * sx, k * (sy = _ ? -1 : 1), tx, ty), reset()) : sy < 0;
    },
    fitExtent: function (extent, object) {
      return (0, _fit.fitExtent)(projection, extent, object);
    },
    fitSize: function (size, object) {
      return (0, _fit.fitSize)(projection, size, object);
    },
    fitWidth: function (width, object) {
      return (0, _fit.fitWidth)(projection, width, object);
    },
    fitHeight: function (height, object) {
      return (0, _fit.fitHeight)(projection, height, object);
    }
  };
}
},{"../clip/rectangle":"../node_modules/d3-geo/src/clip/rectangle.js","../identity":"../node_modules/d3-geo/src/identity.js","../transform":"../node_modules/d3-geo/src/transform.js","./fit":"../node_modules/d3-geo/src/projection/fit.js"}],"../node_modules/d3-geo/src/projection/naturalEarth1.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.naturalEarth1Raw = naturalEarth1Raw;
exports.default = _default;

var _index = _interopRequireDefault(require("./index"));

var _math = require("../math");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function naturalEarth1Raw(lambda, phi) {
  var phi2 = phi * phi,
      phi4 = phi2 * phi2;
  return [lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (0.003971 * phi2 - 0.001529 * phi4))), phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4)))];
}

naturalEarth1Raw.invert = function (x, y) {
  var phi = y,
      i = 25,
      delta;

  do {
    var phi2 = phi * phi,
        phi4 = phi2 * phi2;
    phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4))) - y) / (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 0.005916 * 11 * phi4)));
  } while ((0, _math.abs)(delta) > _math.epsilon && --i > 0);

  return [x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (0.003971 - 0.001529 * phi2)))), phi];
};

function _default() {
  return (0, _index.default)(naturalEarth1Raw).scale(175.295);
}
},{"./index":"../node_modules/d3-geo/src/projection/index.js","../math":"../node_modules/d3-geo/src/math.js"}],"../node_modules/d3-geo/src/projection/orthographic.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.orthographicRaw = orthographicRaw;
exports.default = _default;

var _math = require("../math");

var _azimuthal = require("./azimuthal");

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function orthographicRaw(x, y) {
  return [(0, _math.cos)(y) * (0, _math.sin)(x), (0, _math.sin)(y)];
}

orthographicRaw.invert = (0, _azimuthal.azimuthalInvert)(_math.asin);

function _default() {
  return (0, _index.default)(orthographicRaw).scale(249.5).clipAngle(90 + _math.epsilon);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./azimuthal":"../node_modules/d3-geo/src/projection/azimuthal.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/stereographic.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.stereographicRaw = stereographicRaw;
exports.default = _default;

var _math = require("../math");

var _azimuthal = require("./azimuthal");

var _index = _interopRequireDefault(require("./index"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function stereographicRaw(x, y) {
  var cy = (0, _math.cos)(y),
      k = 1 + (0, _math.cos)(x) * cy;
  return [cy * (0, _math.sin)(x) / k, (0, _math.sin)(y) / k];
}

stereographicRaw.invert = (0, _azimuthal.azimuthalInvert)(function (z) {
  return 2 * (0, _math.atan)(z);
});

function _default() {
  return (0, _index.default)(stereographicRaw).scale(250).clipAngle(142);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./azimuthal":"../node_modules/d3-geo/src/projection/azimuthal.js","./index":"../node_modules/d3-geo/src/projection/index.js"}],"../node_modules/d3-geo/src/projection/transverseMercator.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.transverseMercatorRaw = transverseMercatorRaw;
exports.default = _default;

var _math = require("../math");

var _mercator = require("./mercator");

function transverseMercatorRaw(lambda, phi) {
  return [(0, _math.log)((0, _math.tan)((_math.halfPi + phi) / 2)), -lambda];
}

transverseMercatorRaw.invert = function (x, y) {
  return [-y, 2 * (0, _math.atan)((0, _math.exp)(x)) - _math.halfPi];
};

function _default() {
  var m = (0, _mercator.mercatorProjection)(transverseMercatorRaw),
      center = m.center,
      rotate = m.rotate;

  m.center = function (_) {
    return arguments.length ? center([-_[1], _[0]]) : (_ = center(), [_[1], -_[0]]);
  };

  m.rotate = function (_) {
    return arguments.length ? rotate([_[0], _[1], _.length > 2 ? _[2] + 90 : 90]) : (_ = rotate(), [_[0], _[1], _[2] - 90]);
  };

  return rotate([0, 0, 90]).scale(159.155);
}
},{"../math":"../node_modules/d3-geo/src/math.js","./mercator":"../node_modules/d3-geo/src/projection/mercator.js"}],"../node_modules/d3-geo/src/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
Object.defineProperty(exports, "geoArea", {
  enumerable: true,
  get: function () {
    return _area.default;
  }
});
Object.defineProperty(exports, "geoBounds", {
  enumerable: true,
  get: function () {
    return _bounds.default;
  }
});
Object.defineProperty(exports, "geoCentroid", {
  enumerable: true,
  get: function () {
    return _centroid.default;
  }
});
Object.defineProperty(exports, "geoCircle", {
  enumerable: true,
  get: function () {
    return _circle.default;
  }
});
Object.defineProperty(exports, "geoClipAntimeridian", {
  enumerable: true,
  get: function () {
    return _antimeridian.default;
  }
});
Object.defineProperty(exports, "geoClipCircle", {
  enumerable: true,
  get: function () {
    return _circle2.default;
  }
});
Object.defineProperty(exports, "geoClipExtent", {
  enumerable: true,
  get: function () {
    return _extent.default;
  }
});
Object.defineProperty(exports, "geoClipRectangle", {
  enumerable: true,
  get: function () {
    return _rectangle.default;
  }
});
Object.defineProperty(exports, "geoContains", {
  enumerable: true,
  get: function () {
    return _contains.default;
  }
});
Object.defineProperty(exports, "geoDistance", {
  enumerable: true,
  get: function () {
    return _distance.default;
  }
});
Object.defineProperty(exports, "geoGraticule", {
  enumerable: true,
  get: function () {
    return _graticule.default;
  }
});
Object.defineProperty(exports, "geoGraticule10", {
  enumerable: true,
  get: function () {
    return _graticule.graticule10;
  }
});
Object.defineProperty(exports, "geoInterpolate", {
  enumerable: true,
  get: function () {
    return _interpolate.default;
  }
});
Object.defineProperty(exports, "geoLength", {
  enumerable: true,
  get: function () {
    return _length.default;
  }
});
Object.defineProperty(exports, "geoPath", {
  enumerable: true,
  get: function () {
    return _index.default;
  }
});
Object.defineProperty(exports, "geoAlbers", {
  enumerable: true,
  get: function () {
    return _albers.default;
  }
});
Object.defineProperty(exports, "geoAlbersUsa", {
  enumerable: true,
  get: function () {
    return _albersUsa.default;
  }
});
Object.defineProperty(exports, "geoAzimuthalEqualArea", {
  enumerable: true,
  get: function () {
    return _azimuthalEqualArea.default;
  }
});
Object.defineProperty(exports, "geoAzimuthalEqualAreaRaw", {
  enumerable: true,
  get: function () {
    return _azimuthalEqualArea.azimuthalEqualAreaRaw;
  }
});
Object.defineProperty(exports, "geoAzimuthalEquidistant", {
  enumerable: true,
  get: function () {
    return _azimuthalEquidistant.default;
  }
});
Object.defineProperty(exports, "geoAzimuthalEquidistantRaw", {
  enumerable: true,
  get: function () {
    return _azimuthalEquidistant.azimuthalEquidistantRaw;
  }
});
Object.defineProperty(exports, "geoConicConformal", {
  enumerable: true,
  get: function () {
    return _conicConformal.default;
  }
});
Object.defineProperty(exports, "geoConicConformalRaw", {
  enumerable: true,
  get: function () {
    return _conicConformal.conicConformalRaw;
  }
});
Object.defineProperty(exports, "geoConicEqualArea", {
  enumerable: true,
  get: function () {
    return _conicEqualArea.default;
  }
});
Object.defineProperty(exports, "geoConicEqualAreaRaw", {
  enumerable: true,
  get: function () {
    return _conicEqualArea.conicEqualAreaRaw;
  }
});
Object.defineProperty(exports, "geoConicEquidistant", {
  enumerable: true,
  get: function () {
    return _conicEquidistant.default;
  }
});
Object.defineProperty(exports, "geoConicEquidistantRaw", {
  enumerable: true,
  get: function () {
    return _conicEquidistant.conicEquidistantRaw;
  }
});
Object.defineProperty(exports, "geoEqualEarth", {
  enumerable: true,
  get: function () {
    return _equalEarth.default;
  }
});
Object.defineProperty(exports, "geoEqualEarthRaw", {
  enumerable: true,
  get: function () {
    return _equalEarth.equalEarthRaw;
  }
});
Object.defineProperty(exports, "geoEquirectangular", {
  enumerable: true,
  get: function () {
    return _equirectangular.default;
  }
});
Object.defineProperty(exports, "geoEquirectangularRaw", {
  enumerable: true,
  get: function () {
    return _equirectangular.equirectangularRaw;
  }
});
Object.defineProperty(exports, "geoGnomonic", {
  enumerable: true,
  get: function () {
    return _gnomonic.default;
  }
});
Object.defineProperty(exports, "geoGnomonicRaw", {
  enumerable: true,
  get: function () {
    return _gnomonic.gnomonicRaw;
  }
});
Object.defineProperty(exports, "geoIdentity", {
  enumerable: true,
  get: function () {
    return _identity.default;
  }
});
Object.defineProperty(exports, "geoProjection", {
  enumerable: true,
  get: function () {
    return _index2.default;
  }
});
Object.defineProperty(exports, "geoProjectionMutator", {
  enumerable: true,
  get: function () {
    return _index2.projectionMutator;
  }
});
Object.defineProperty(exports, "geoMercator", {
  enumerable: true,
  get: function () {
    return _mercator.default;
  }
});
Object.defineProperty(exports, "geoMercatorRaw", {
  enumerable: true,
  get: function () {
    return _mercator.mercatorRaw;
  }
});
Object.defineProperty(exports, "geoNaturalEarth1", {
  enumerable: true,
  get: function () {
    return _naturalEarth.default;
  }
});
Object.defineProperty(exports, "geoNaturalEarth1Raw", {
  enumerable: true,
  get: function () {
    return _naturalEarth.naturalEarth1Raw;
  }
});
Object.defineProperty(exports, "geoOrthographic", {
  enumerable: true,
  get: function () {
    return _orthographic.default;
  }
});
Object.defineProperty(exports, "geoOrthographicRaw", {
  enumerable: true,
  get: function () {
    return _orthographic.orthographicRaw;
  }
});
Object.defineProperty(exports, "geoStereographic", {
  enumerable: true,
  get: function () {
    return _stereographic.default;
  }
});
Object.defineProperty(exports, "geoStereographicRaw", {
  enumerable: true,
  get: function () {
    return _stereographic.stereographicRaw;
  }
});
Object.defineProperty(exports, "geoTransverseMercator", {
  enumerable: true,
  get: function () {
    return _transverseMercator.default;
  }
});
Object.defineProperty(exports, "geoTransverseMercatorRaw", {
  enumerable: true,
  get: function () {
    return _transverseMercator.transverseMercatorRaw;
  }
});
Object.defineProperty(exports, "geoRotation", {
  enumerable: true,
  get: function () {
    return _rotation.default;
  }
});
Object.defineProperty(exports, "geoStream", {
  enumerable: true,
  get: function () {
    return _stream.default;
  }
});
Object.defineProperty(exports, "geoTransform", {
  enumerable: true,
  get: function () {
    return _transform.default;
  }
});

var _area = _interopRequireDefault(require("./area"));

var _bounds = _interopRequireDefault(require("./bounds"));

var _centroid = _interopRequireDefault(require("./centroid"));

var _circle = _interopRequireDefault(require("./circle"));

var _antimeridian = _interopRequireDefault(require("./clip/antimeridian"));

var _circle2 = _interopRequireDefault(require("./clip/circle"));

var _extent = _interopRequireDefault(require("./clip/extent"));

var _rectangle = _interopRequireDefault(require("./clip/rectangle"));

var _contains = _interopRequireDefault(require("./contains"));

var _distance = _interopRequireDefault(require("./distance"));

var _graticule = _interopRequireWildcard(require("./graticule"));

var _interpolate = _interopRequireDefault(require("./interpolate"));

var _length = _interopRequireDefault(require("./length"));

var _index = _interopRequireDefault(require("./path/index"));

var _albers = _interopRequireDefault(require("./projection/albers"));

var _albersUsa = _interopRequireDefault(require("./projection/albersUsa"));

var _azimuthalEqualArea = _interopRequireWildcard(require("./projection/azimuthalEqualArea"));

var _azimuthalEquidistant = _interopRequireWildcard(require("./projection/azimuthalEquidistant"));

var _conicConformal = _interopRequireWildcard(require("./projection/conicConformal"));

var _conicEqualArea = _interopRequireWildcard(require("./projection/conicEqualArea"));

var _conicEquidistant = _interopRequireWildcard(require("./projection/conicEquidistant"));

var _equalEarth = _interopRequireWildcard(require("./projection/equalEarth"));

var _equirectangular = _interopRequireWildcard(require("./projection/equirectangular"));

var _gnomonic = _interopRequireWildcard(require("./projection/gnomonic"));

var _identity = _interopRequireDefault(require("./projection/identity"));

var _index2 = _interopRequireWildcard(require("./projection/index"));

var _mercator = _interopRequireWildcard(require("./projection/mercator"));

var _naturalEarth = _interopRequireWildcard(require("./projection/naturalEarth1"));

var _orthographic = _interopRequireWildcard(require("./projection/orthographic"));

var _stereographic = _interopRequireWildcard(require("./projection/stereographic"));

var _transverseMercator = _interopRequireWildcard(require("./projection/transverseMercator"));

var _rotation = _interopRequireDefault(require("./rotation"));

var _stream = _interopRequireDefault(require("./stream"));

var _transform = _interopRequireDefault(require("./transform"));

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = Object.defineProperty && Object.getOwnPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : {}; if (desc.get || desc.set) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }
},{"./area":"../node_modules/d3-geo/src/area.js","./bounds":"../node_modules/d3-geo/src/bounds.js","./centroid":"../node_modules/d3-geo/src/centroid.js","./circle":"../node_modules/d3-geo/src/circle.js","./clip/antimeridian":"../node_modules/d3-geo/src/clip/antimeridian.js","./clip/circle":"../node_modules/d3-geo/src/clip/circle.js","./clip/extent":"../node_modules/d3-geo/src/clip/extent.js","./clip/rectangle":"../node_modules/d3-geo/src/clip/rectangle.js","./contains":"../node_modules/d3-geo/src/contains.js","./distance":"../node_modules/d3-geo/src/distance.js","./graticule":"../node_modules/d3-geo/src/graticule.js","./interpolate":"../node_modules/d3-geo/src/interpolate.js","./length":"../node_modules/d3-geo/src/length.js","./path/index":"../node_modules/d3-geo/src/path/index.js","./projection/albers":"../node_modules/d3-geo/src/projection/albers.js","./projection/albersUsa":"../node_modules/d3-geo/src/projection/albersUsa.js","./projection/azimuthalEqualArea":"../node_modules/d3-geo/src/projection/azimuthalEqualArea.js","./projection/azimuthalEquidistant":"../node_modules/d3-geo/src/projection/azimuthalEquidistant.js","./projection/conicConformal":"../node_modules/d3-geo/src/projection/conicConformal.js","./projection/conicEqualArea":"../node_modules/d3-geo/src/projection/conicEqualArea.js","./projection/conicEquidistant":"../node_modules/d3-geo/src/projection/conicEquidistant.js","./projection/equalEarth":"../node_modules/d3-geo/src/projection/equalEarth.js","./projection/equirectangular":"../node_modules/d3-geo/src/projection/equirectangular.js","./projection/gnomonic":"../node_modules/d3-geo/src/projection/gnomonic.js","./projection/identity":"../node_modules/d3-geo/src/projection/identity.js","./projection/index":"../node_modules/d3-geo/src/projection/index.js","./projection/mercator":"../node_modules/d3-geo/src/projection/mercator.js","./projection/naturalEarth1":"../node_modules/d3-geo/src/projection/naturalEarth1.js","./projection/orthographic":"../node_modules/d3-geo/src/projection/orthographic.js","./projection/stereographic":"../node_modules/d3-geo/src/projection/stereographic.js","./projection/transverseMercator":"../node_modules/d3-geo/src/projection/transverseMercator.js","./rotation":"../node_modules/d3-geo/src/rotation.js","./stream":"../node_modules/d3-geo/src/stream.js","./transform":"../node_modules/d3-geo/src/transform.js"}],"../node_modules/d3-selection/src/namespaces.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = exports.xhtml = void 0;
var xhtml = "http://www.w3.org/1999/xhtml";
exports.xhtml = xhtml;
var _default = {
  svg: "http://www.w3.org/2000/svg",
  xhtml: xhtml,
  xlink: "http://www.w3.org/1999/xlink",
  xml: "http://www.w3.org/XML/1998/namespace",
  xmlns: "http://www.w3.org/2000/xmlns/"
};
exports.default = _default;
},{}],"../node_modules/d3-selection/src/namespace.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _namespaces = _interopRequireDefault(require("./namespaces"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(name) {
  var prefix = name += "",
      i = prefix.indexOf(":");
  if (i >= 0 && (prefix = name.slice(0, i)) !== "xmlns") name = name.slice(i + 1);
  return _namespaces.default.hasOwnProperty(prefix) ? {
    space: _namespaces.default[prefix],
    local: name
  } : name;
}
},{"./namespaces":"../node_modules/d3-selection/src/namespaces.js"}],"../node_modules/d3-selection/src/creator.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _namespace = _interopRequireDefault(require("./namespace"));

var _namespaces = require("./namespaces");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function creatorInherit(name) {
  return function () {
    var document = this.ownerDocument,
        uri = this.namespaceURI;
    return uri === _namespaces.xhtml && document.documentElement.namespaceURI === _namespaces.xhtml ? document.createElement(name) : document.createElementNS(uri, name);
  };
}

function creatorFixed(fullname) {
  return function () {
    return this.ownerDocument.createElementNS(fullname.space, fullname.local);
  };
}

function _default(name) {
  var fullname = (0, _namespace.default)(name);
  return (fullname.local ? creatorFixed : creatorInherit)(fullname);
}
},{"./namespace":"../node_modules/d3-selection/src/namespace.js","./namespaces":"../node_modules/d3-selection/src/namespaces.js"}],"../node_modules/d3-selection/src/selector.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function none() {}

function _default(selector) {
  return selector == null ? none : function () {
    return this.querySelector(selector);
  };
}
},{}],"../node_modules/d3-selection/src/selection/select.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

var _selector = _interopRequireDefault(require("../selector"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(select) {
  if (typeof select !== "function") select = (0, _selector.default)(select);

  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = new Array(n), node, subnode, i = 0; i < n; ++i) {
      if ((node = group[i]) && (subnode = select.call(node, node.__data__, i, group))) {
        if ("__data__" in node) subnode.__data__ = node.__data__;
        subgroup[i] = subnode;
      }
    }
  }

  return new _index.Selection(subgroups, this._parents);
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js","../selector":"../node_modules/d3-selection/src/selector.js"}],"../node_modules/d3-selection/src/selectorAll.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function empty() {
  return [];
}

function _default(selector) {
  return selector == null ? empty : function () {
    return this.querySelectorAll(selector);
  };
}
},{}],"../node_modules/d3-selection/src/selection/selectAll.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

var _selectorAll = _interopRequireDefault(require("../selectorAll"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(select) {
  if (typeof select !== "function") select = (0, _selectorAll.default)(select);

  for (var groups = this._groups, m = groups.length, subgroups = [], parents = [], j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        subgroups.push(select.call(node, node.__data__, i, group));
        parents.push(node);
      }
    }
  }

  return new _index.Selection(subgroups, parents);
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js","../selectorAll":"../node_modules/d3-selection/src/selectorAll.js"}],"../node_modules/d3-selection/src/matcher.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(selector) {
  return function () {
    return this.matches(selector);
  };
}
},{}],"../node_modules/d3-selection/src/selection/filter.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

var _matcher = _interopRequireDefault(require("../matcher"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(match) {
  if (typeof match !== "function") match = (0, _matcher.default)(match);

  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = [], node, i = 0; i < n; ++i) {
      if ((node = group[i]) && match.call(node, node.__data__, i, group)) {
        subgroup.push(node);
      }
    }
  }

  return new _index.Selection(subgroups, this._parents);
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js","../matcher":"../node_modules/d3-selection/src/matcher.js"}],"../node_modules/d3-selection/src/selection/sparse.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(update) {
  return new Array(update.length);
}
},{}],"../node_modules/d3-selection/src/selection/enter.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.EnterNode = EnterNode;

var _sparse = _interopRequireDefault(require("./sparse"));

var _index = require("./index");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  return new _index.Selection(this._enter || this._groups.map(_sparse.default), this._parents);
}

function EnterNode(parent, datum) {
  this.ownerDocument = parent.ownerDocument;
  this.namespaceURI = parent.namespaceURI;
  this._next = null;
  this._parent = parent;
  this.__data__ = datum;
}

EnterNode.prototype = {
  constructor: EnterNode,
  appendChild: function (child) {
    return this._parent.insertBefore(child, this._next);
  },
  insertBefore: function (child, next) {
    return this._parent.insertBefore(child, next);
  },
  querySelector: function (selector) {
    return this._parent.querySelector(selector);
  },
  querySelectorAll: function (selector) {
    return this._parent.querySelectorAll(selector);
  }
};
},{"./sparse":"../node_modules/d3-selection/src/selection/sparse.js","./index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/constant.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(x) {
  return function () {
    return x;
  };
}
},{}],"../node_modules/d3-selection/src/selection/data.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

var _enter = require("./enter");

var _constant = _interopRequireDefault(require("../constant"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var keyPrefix = "$"; // Protect against keys like “__proto__”.

function bindIndex(parent, group, enter, update, exit, data) {
  var i = 0,
      node,
      groupLength = group.length,
      dataLength = data.length; // Put any non-null nodes that fit into update.
  // Put any null nodes into enter.
  // Put any remaining data into enter.

  for (; i < dataLength; ++i) {
    if (node = group[i]) {
      node.__data__ = data[i];
      update[i] = node;
    } else {
      enter[i] = new _enter.EnterNode(parent, data[i]);
    }
  } // Put any non-null nodes that don’t fit into exit.


  for (; i < groupLength; ++i) {
    if (node = group[i]) {
      exit[i] = node;
    }
  }
}

function bindKey(parent, group, enter, update, exit, data, key) {
  var i,
      node,
      nodeByKeyValue = {},
      groupLength = group.length,
      dataLength = data.length,
      keyValues = new Array(groupLength),
      keyValue; // Compute the key for each node.
  // If multiple nodes have the same key, the duplicates are added to exit.

  for (i = 0; i < groupLength; ++i) {
    if (node = group[i]) {
      keyValues[i] = keyValue = keyPrefix + key.call(node, node.__data__, i, group);

      if (keyValue in nodeByKeyValue) {
        exit[i] = node;
      } else {
        nodeByKeyValue[keyValue] = node;
      }
    }
  } // Compute the key for each datum.
  // If there a node associated with this key, join and add it to update.
  // If there is not (or the key is a duplicate), add it to enter.


  for (i = 0; i < dataLength; ++i) {
    keyValue = keyPrefix + key.call(parent, data[i], i, data);

    if (node = nodeByKeyValue[keyValue]) {
      update[i] = node;
      node.__data__ = data[i];
      nodeByKeyValue[keyValue] = null;
    } else {
      enter[i] = new _enter.EnterNode(parent, data[i]);
    }
  } // Add any remaining nodes that were not bound to data to exit.


  for (i = 0; i < groupLength; ++i) {
    if ((node = group[i]) && nodeByKeyValue[keyValues[i]] === node) {
      exit[i] = node;
    }
  }
}

function _default(value, key) {
  if (!value) {
    data = new Array(this.size()), j = -1;
    this.each(function (d) {
      data[++j] = d;
    });
    return data;
  }

  var bind = key ? bindKey : bindIndex,
      parents = this._parents,
      groups = this._groups;
  if (typeof value !== "function") value = (0, _constant.default)(value);

  for (var m = groups.length, update = new Array(m), enter = new Array(m), exit = new Array(m), j = 0; j < m; ++j) {
    var parent = parents[j],
        group = groups[j],
        groupLength = group.length,
        data = value.call(parent, parent && parent.__data__, j, parents),
        dataLength = data.length,
        enterGroup = enter[j] = new Array(dataLength),
        updateGroup = update[j] = new Array(dataLength),
        exitGroup = exit[j] = new Array(groupLength);
    bind(parent, group, enterGroup, updateGroup, exitGroup, data, key); // Now connect the enter nodes to their following update node, such that
    // appendChild can insert the materialized enter node before this node,
    // rather than at the end of the parent node.

    for (var i0 = 0, i1 = 0, previous, next; i0 < dataLength; ++i0) {
      if (previous = enterGroup[i0]) {
        if (i0 >= i1) i1 = i0 + 1;

        while (!(next = updateGroup[i1]) && ++i1 < dataLength);

        previous._next = next || null;
      }
    }
  }

  update = new _index.Selection(update, parents);
  update._enter = enter;
  update._exit = exit;
  return update;
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js","./enter":"../node_modules/d3-selection/src/selection/enter.js","../constant":"../node_modules/d3-selection/src/constant.js"}],"../node_modules/d3-selection/src/selection/exit.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _sparse = _interopRequireDefault(require("./sparse"));

var _index = require("./index");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default() {
  return new _index.Selection(this._exit || this._groups.map(_sparse.default), this._parents);
}
},{"./sparse":"../node_modules/d3-selection/src/selection/sparse.js","./index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/selection/join.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(onenter, onupdate, onexit) {
  var enter = this.enter(),
      update = this,
      exit = this.exit();
  enter = typeof onenter === "function" ? onenter(enter) : enter.append(onenter + "");
  if (onupdate != null) update = onupdate(update);
  if (onexit == null) exit.remove();else onexit(exit);
  return enter && update ? enter.merge(update).order() : update;
}
},{}],"../node_modules/d3-selection/src/selection/merge.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

function _default(selection) {
  for (var groups0 = this._groups, groups1 = selection._groups, m0 = groups0.length, m1 = groups1.length, m = Math.min(m0, m1), merges = new Array(m0), j = 0; j < m; ++j) {
    for (var group0 = groups0[j], group1 = groups1[j], n = group0.length, merge = merges[j] = new Array(n), node, i = 0; i < n; ++i) {
      if (node = group0[i] || group1[i]) {
        merge[i] = node;
      }
    }
  }

  for (; j < m0; ++j) {
    merges[j] = groups0[j];
  }

  return new _index.Selection(merges, this._parents);
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/selection/order.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  for (var groups = this._groups, j = -1, m = groups.length; ++j < m;) {
    for (var group = groups[j], i = group.length - 1, next = group[i], node; --i >= 0;) {
      if (node = group[i]) {
        if (next && node.compareDocumentPosition(next) ^ 4) next.parentNode.insertBefore(node, next);
        next = node;
      }
    }
  }

  return this;
}
},{}],"../node_modules/d3-selection/src/selection/sort.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./index");

function _default(compare) {
  if (!compare) compare = ascending;

  function compareNode(a, b) {
    return a && b ? compare(a.__data__, b.__data__) : !a - !b;
  }

  for (var groups = this._groups, m = groups.length, sortgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, sortgroup = sortgroups[j] = new Array(n), node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        sortgroup[i] = node;
      }
    }

    sortgroup.sort(compareNode);
  }

  return new _index.Selection(sortgroups, this._parents).order();
}

function ascending(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}
},{"./index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/selection/call.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  var callback = arguments[0];
  arguments[0] = this;
  callback.apply(null, arguments);
  return this;
}
},{}],"../node_modules/d3-selection/src/selection/nodes.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  var nodes = new Array(this.size()),
      i = -1;
  this.each(function () {
    nodes[++i] = this;
  });
  return nodes;
}
},{}],"../node_modules/d3-selection/src/selection/node.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
    for (var group = groups[j], i = 0, n = group.length; i < n; ++i) {
      var node = group[i];
      if (node) return node;
    }
  }

  return null;
}
},{}],"../node_modules/d3-selection/src/selection/size.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  var size = 0;
  this.each(function () {
    ++size;
  });
  return size;
}
},{}],"../node_modules/d3-selection/src/selection/empty.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default() {
  return !this.node();
}
},{}],"../node_modules/d3-selection/src/selection/each.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(callback) {
  for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
    for (var group = groups[j], i = 0, n = group.length, node; i < n; ++i) {
      if (node = group[i]) callback.call(node, node.__data__, i, group);
    }
  }

  return this;
}
},{}],"../node_modules/d3-selection/src/selection/attr.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _namespace = _interopRequireDefault(require("../namespace"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function attrRemove(name) {
  return function () {
    this.removeAttribute(name);
  };
}

function attrRemoveNS(fullname) {
  return function () {
    this.removeAttributeNS(fullname.space, fullname.local);
  };
}

function attrConstant(name, value) {
  return function () {
    this.setAttribute(name, value);
  };
}

function attrConstantNS(fullname, value) {
  return function () {
    this.setAttributeNS(fullname.space, fullname.local, value);
  };
}

function attrFunction(name, value) {
  return function () {
    var v = value.apply(this, arguments);
    if (v == null) this.removeAttribute(name);else this.setAttribute(name, v);
  };
}

function attrFunctionNS(fullname, value) {
  return function () {
    var v = value.apply(this, arguments);
    if (v == null) this.removeAttributeNS(fullname.space, fullname.local);else this.setAttributeNS(fullname.space, fullname.local, v);
  };
}

function _default(name, value) {
  var fullname = (0, _namespace.default)(name);

  if (arguments.length < 2) {
    var node = this.node();
    return fullname.local ? node.getAttributeNS(fullname.space, fullname.local) : node.getAttribute(fullname);
  }

  return this.each((value == null ? fullname.local ? attrRemoveNS : attrRemove : typeof value === "function" ? fullname.local ? attrFunctionNS : attrFunction : fullname.local ? attrConstantNS : attrConstant)(fullname, value));
}
},{"../namespace":"../node_modules/d3-selection/src/namespace.js"}],"../node_modules/d3-selection/src/window.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(node) {
  return node.ownerDocument && node.ownerDocument.defaultView || // node is a Node
  node.document && node // node is a Window
  || node.defaultView; // node is a Document
}
},{}],"../node_modules/d3-selection/src/selection/style.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.styleValue = styleValue;

var _window = _interopRequireDefault(require("../window"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function styleRemove(name) {
  return function () {
    this.style.removeProperty(name);
  };
}

function styleConstant(name, value, priority) {
  return function () {
    this.style.setProperty(name, value, priority);
  };
}

function styleFunction(name, value, priority) {
  return function () {
    var v = value.apply(this, arguments);
    if (v == null) this.style.removeProperty(name);else this.style.setProperty(name, v, priority);
  };
}

function _default(name, value, priority) {
  return arguments.length > 1 ? this.each((value == null ? styleRemove : typeof value === "function" ? styleFunction : styleConstant)(name, value, priority == null ? "" : priority)) : styleValue(this.node(), name);
}

function styleValue(node, name) {
  return node.style.getPropertyValue(name) || (0, _window.default)(node).getComputedStyle(node, null).getPropertyValue(name);
}
},{"../window":"../node_modules/d3-selection/src/window.js"}],"../node_modules/d3-selection/src/selection/property.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function propertyRemove(name) {
  return function () {
    delete this[name];
  };
}

function propertyConstant(name, value) {
  return function () {
    this[name] = value;
  };
}

function propertyFunction(name, value) {
  return function () {
    var v = value.apply(this, arguments);
    if (v == null) delete this[name];else this[name] = v;
  };
}

function _default(name, value) {
  return arguments.length > 1 ? this.each((value == null ? propertyRemove : typeof value === "function" ? propertyFunction : propertyConstant)(name, value)) : this.node()[name];
}
},{}],"../node_modules/d3-selection/src/selection/classed.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function classArray(string) {
  return string.trim().split(/^|\s+/);
}

function classList(node) {
  return node.classList || new ClassList(node);
}

function ClassList(node) {
  this._node = node;
  this._names = classArray(node.getAttribute("class") || "");
}

ClassList.prototype = {
  add: function (name) {
    var i = this._names.indexOf(name);

    if (i < 0) {
      this._names.push(name);

      this._node.setAttribute("class", this._names.join(" "));
    }
  },
  remove: function (name) {
    var i = this._names.indexOf(name);

    if (i >= 0) {
      this._names.splice(i, 1);

      this._node.setAttribute("class", this._names.join(" "));
    }
  },
  contains: function (name) {
    return this._names.indexOf(name) >= 0;
  }
};

function classedAdd(node, names) {
  var list = classList(node),
      i = -1,
      n = names.length;

  while (++i < n) list.add(names[i]);
}

function classedRemove(node, names) {
  var list = classList(node),
      i = -1,
      n = names.length;

  while (++i < n) list.remove(names[i]);
}

function classedTrue(names) {
  return function () {
    classedAdd(this, names);
  };
}

function classedFalse(names) {
  return function () {
    classedRemove(this, names);
  };
}

function classedFunction(names, value) {
  return function () {
    (value.apply(this, arguments) ? classedAdd : classedRemove)(this, names);
  };
}

function _default(name, value) {
  var names = classArray(name + "");

  if (arguments.length < 2) {
    var list = classList(this.node()),
        i = -1,
        n = names.length;

    while (++i < n) if (!list.contains(names[i])) return false;

    return true;
  }

  return this.each((typeof value === "function" ? classedFunction : value ? classedTrue : classedFalse)(names, value));
}
},{}],"../node_modules/d3-selection/src/selection/text.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function textRemove() {
  this.textContent = "";
}

function textConstant(value) {
  return function () {
    this.textContent = value;
  };
}

function textFunction(value) {
  return function () {
    var v = value.apply(this, arguments);
    this.textContent = v == null ? "" : v;
  };
}

function _default(value) {
  return arguments.length ? this.each(value == null ? textRemove : (typeof value === "function" ? textFunction : textConstant)(value)) : this.node().textContent;
}
},{}],"../node_modules/d3-selection/src/selection/html.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function htmlRemove() {
  this.innerHTML = "";
}

function htmlConstant(value) {
  return function () {
    this.innerHTML = value;
  };
}

function htmlFunction(value) {
  return function () {
    var v = value.apply(this, arguments);
    this.innerHTML = v == null ? "" : v;
  };
}

function _default(value) {
  return arguments.length ? this.each(value == null ? htmlRemove : (typeof value === "function" ? htmlFunction : htmlConstant)(value)) : this.node().innerHTML;
}
},{}],"../node_modules/d3-selection/src/selection/raise.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function raise() {
  if (this.nextSibling) this.parentNode.appendChild(this);
}

function _default() {
  return this.each(raise);
}
},{}],"../node_modules/d3-selection/src/selection/lower.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function lower() {
  if (this.previousSibling) this.parentNode.insertBefore(this, this.parentNode.firstChild);
}

function _default() {
  return this.each(lower);
}
},{}],"../node_modules/d3-selection/src/selection/append.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _creator = _interopRequireDefault(require("../creator"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(name) {
  var create = typeof name === "function" ? name : (0, _creator.default)(name);
  return this.select(function () {
    return this.appendChild(create.apply(this, arguments));
  });
}
},{"../creator":"../node_modules/d3-selection/src/creator.js"}],"../node_modules/d3-selection/src/selection/insert.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _creator = _interopRequireDefault(require("../creator"));

var _selector = _interopRequireDefault(require("../selector"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function constantNull() {
  return null;
}

function _default(name, before) {
  var create = typeof name === "function" ? name : (0, _creator.default)(name),
      select = before == null ? constantNull : typeof before === "function" ? before : (0, _selector.default)(before);
  return this.select(function () {
    return this.insertBefore(create.apply(this, arguments), select.apply(this, arguments) || null);
  });
}
},{"../creator":"../node_modules/d3-selection/src/creator.js","../selector":"../node_modules/d3-selection/src/selector.js"}],"../node_modules/d3-selection/src/selection/remove.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function remove() {
  var parent = this.parentNode;
  if (parent) parent.removeChild(this);
}

function _default() {
  return this.each(remove);
}
},{}],"../node_modules/d3-selection/src/selection/clone.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function selection_cloneShallow() {
  return this.parentNode.insertBefore(this.cloneNode(false), this.nextSibling);
}

function selection_cloneDeep() {
  return this.parentNode.insertBefore(this.cloneNode(true), this.nextSibling);
}

function _default(deep) {
  return this.select(deep ? selection_cloneDeep : selection_cloneShallow);
}
},{}],"../node_modules/d3-selection/src/selection/datum.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(value) {
  return arguments.length ? this.property("__data__", value) : this.node().__data__;
}
},{}],"../node_modules/d3-selection/src/selection/on.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;
exports.customEvent = customEvent;
exports.event = void 0;
var filterEvents = {};
var event = null;
exports.event = event;

if (typeof document !== "undefined") {
  var element = document.documentElement;

  if (!("onmouseenter" in element)) {
    filterEvents = {
      mouseenter: "mouseover",
      mouseleave: "mouseout"
    };
  }
}

function filterContextListener(listener, index, group) {
  listener = contextListener(listener, index, group);
  return function (event) {
    var related = event.relatedTarget;

    if (!related || related !== this && !(related.compareDocumentPosition(this) & 8)) {
      listener.call(this, event);
    }
  };
}

function contextListener(listener, index, group) {
  return function (event1) {
    var event0 = event; // Events can be reentrant (e.g., focus).

    exports.event = event = event1;

    try {
      listener.call(this, this.__data__, index, group);
    } finally {
      exports.event = event = event0;
    }
  };
}

function parseTypenames(typenames) {
  return typenames.trim().split(/^|\s+/).map(function (t) {
    var name = "",
        i = t.indexOf(".");
    if (i >= 0) name = t.slice(i + 1), t = t.slice(0, i);
    return {
      type: t,
      name: name
    };
  });
}

function onRemove(typename) {
  return function () {
    var on = this.__on;
    if (!on) return;

    for (var j = 0, i = -1, m = on.length, o; j < m; ++j) {
      if (o = on[j], (!typename.type || o.type === typename.type) && o.name === typename.name) {
        this.removeEventListener(o.type, o.listener, o.capture);
      } else {
        on[++i] = o;
      }
    }

    if (++i) on.length = i;else delete this.__on;
  };
}

function onAdd(typename, value, capture) {
  var wrap = filterEvents.hasOwnProperty(typename.type) ? filterContextListener : contextListener;
  return function (d, i, group) {
    var on = this.__on,
        o,
        listener = wrap(value, i, group);
    if (on) for (var j = 0, m = on.length; j < m; ++j) {
      if ((o = on[j]).type === typename.type && o.name === typename.name) {
        this.removeEventListener(o.type, o.listener, o.capture);
        this.addEventListener(o.type, o.listener = listener, o.capture = capture);
        o.value = value;
        return;
      }
    }
    this.addEventListener(typename.type, listener, capture);
    o = {
      type: typename.type,
      name: typename.name,
      value: value,
      listener: listener,
      capture: capture
    };
    if (!on) this.__on = [o];else on.push(o);
  };
}

function _default(typename, value, capture) {
  var typenames = parseTypenames(typename + ""),
      i,
      n = typenames.length,
      t;

  if (arguments.length < 2) {
    var on = this.node().__on;

    if (on) for (var j = 0, m = on.length, o; j < m; ++j) {
      for (i = 0, o = on[j]; i < n; ++i) {
        if ((t = typenames[i]).type === o.type && t.name === o.name) {
          return o.value;
        }
      }
    }
    return;
  }

  on = value ? onAdd : onRemove;
  if (capture == null) capture = false;

  for (i = 0; i < n; ++i) this.each(on(typenames[i], value, capture));

  return this;
}

function customEvent(event1, listener, that, args) {
  var event0 = event;
  event1.sourceEvent = event;
  exports.event = event = event1;

  try {
    return listener.apply(that, args);
  } finally {
    exports.event = event = event0;
  }
}
},{}],"../node_modules/d3-selection/src/selection/dispatch.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _window = _interopRequireDefault(require("../window"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function dispatchEvent(node, type, params) {
  var window = (0, _window.default)(node),
      event = window.CustomEvent;

  if (typeof event === "function") {
    event = new event(type, params);
  } else {
    event = window.document.createEvent("Event");
    if (params) event.initEvent(type, params.bubbles, params.cancelable), event.detail = params.detail;else event.initEvent(type, false, false);
  }

  node.dispatchEvent(event);
}

function dispatchConstant(type, params) {
  return function () {
    return dispatchEvent(this, type, params);
  };
}

function dispatchFunction(type, params) {
  return function () {
    return dispatchEvent(this, type, params.apply(this, arguments));
  };
}

function _default(type, params) {
  return this.each((typeof params === "function" ? dispatchFunction : dispatchConstant)(type, params));
}
},{"../window":"../node_modules/d3-selection/src/window.js"}],"../node_modules/d3-selection/src/selection/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.Selection = Selection;
exports.default = exports.root = void 0;

var _select = _interopRequireDefault(require("./select"));

var _selectAll = _interopRequireDefault(require("./selectAll"));

var _filter = _interopRequireDefault(require("./filter"));

var _data = _interopRequireDefault(require("./data"));

var _enter = _interopRequireDefault(require("./enter"));

var _exit = _interopRequireDefault(require("./exit"));

var _join = _interopRequireDefault(require("./join"));

var _merge = _interopRequireDefault(require("./merge"));

var _order = _interopRequireDefault(require("./order"));

var _sort = _interopRequireDefault(require("./sort"));

var _call = _interopRequireDefault(require("./call"));

var _nodes = _interopRequireDefault(require("./nodes"));

var _node = _interopRequireDefault(require("./node"));

var _size = _interopRequireDefault(require("./size"));

var _empty = _interopRequireDefault(require("./empty"));

var _each = _interopRequireDefault(require("./each"));

var _attr = _interopRequireDefault(require("./attr"));

var _style = _interopRequireDefault(require("./style"));

var _property = _interopRequireDefault(require("./property"));

var _classed = _interopRequireDefault(require("./classed"));

var _text = _interopRequireDefault(require("./text"));

var _html = _interopRequireDefault(require("./html"));

var _raise = _interopRequireDefault(require("./raise"));

var _lower = _interopRequireDefault(require("./lower"));

var _append = _interopRequireDefault(require("./append"));

var _insert = _interopRequireDefault(require("./insert"));

var _remove = _interopRequireDefault(require("./remove"));

var _clone = _interopRequireDefault(require("./clone"));

var _datum = _interopRequireDefault(require("./datum"));

var _on = _interopRequireDefault(require("./on"));

var _dispatch = _interopRequireDefault(require("./dispatch"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var root = [null];
exports.root = root;

function Selection(groups, parents) {
  this._groups = groups;
  this._parents = parents;
}

function selection() {
  return new Selection([[document.documentElement]], root);
}

Selection.prototype = selection.prototype = {
  constructor: Selection,
  select: _select.default,
  selectAll: _selectAll.default,
  filter: _filter.default,
  data: _data.default,
  enter: _enter.default,
  exit: _exit.default,
  join: _join.default,
  merge: _merge.default,
  order: _order.default,
  sort: _sort.default,
  call: _call.default,
  nodes: _nodes.default,
  node: _node.default,
  size: _size.default,
  empty: _empty.default,
  each: _each.default,
  attr: _attr.default,
  style: _style.default,
  property: _property.default,
  classed: _classed.default,
  text: _text.default,
  html: _html.default,
  raise: _raise.default,
  lower: _lower.default,
  append: _append.default,
  insert: _insert.default,
  remove: _remove.default,
  clone: _clone.default,
  datum: _datum.default,
  on: _on.default,
  dispatch: _dispatch.default
};
var _default = selection;
exports.default = _default;
},{"./select":"../node_modules/d3-selection/src/selection/select.js","./selectAll":"../node_modules/d3-selection/src/selection/selectAll.js","./filter":"../node_modules/d3-selection/src/selection/filter.js","./data":"../node_modules/d3-selection/src/selection/data.js","./enter":"../node_modules/d3-selection/src/selection/enter.js","./exit":"../node_modules/d3-selection/src/selection/exit.js","./join":"../node_modules/d3-selection/src/selection/join.js","./merge":"../node_modules/d3-selection/src/selection/merge.js","./order":"../node_modules/d3-selection/src/selection/order.js","./sort":"../node_modules/d3-selection/src/selection/sort.js","./call":"../node_modules/d3-selection/src/selection/call.js","./nodes":"../node_modules/d3-selection/src/selection/nodes.js","./node":"../node_modules/d3-selection/src/selection/node.js","./size":"../node_modules/d3-selection/src/selection/size.js","./empty":"../node_modules/d3-selection/src/selection/empty.js","./each":"../node_modules/d3-selection/src/selection/each.js","./attr":"../node_modules/d3-selection/src/selection/attr.js","./style":"../node_modules/d3-selection/src/selection/style.js","./property":"../node_modules/d3-selection/src/selection/property.js","./classed":"../node_modules/d3-selection/src/selection/classed.js","./text":"../node_modules/d3-selection/src/selection/text.js","./html":"../node_modules/d3-selection/src/selection/html.js","./raise":"../node_modules/d3-selection/src/selection/raise.js","./lower":"../node_modules/d3-selection/src/selection/lower.js","./append":"../node_modules/d3-selection/src/selection/append.js","./insert":"../node_modules/d3-selection/src/selection/insert.js","./remove":"../node_modules/d3-selection/src/selection/remove.js","./clone":"../node_modules/d3-selection/src/selection/clone.js","./datum":"../node_modules/d3-selection/src/selection/datum.js","./on":"../node_modules/d3-selection/src/selection/on.js","./dispatch":"../node_modules/d3-selection/src/selection/dispatch.js"}],"../node_modules/d3-selection/src/select.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./selection/index");

function _default(selector) {
  return typeof selector === "string" ? new _index.Selection([[document.querySelector(selector)]], [document.documentElement]) : new _index.Selection([[selector]], _index.root);
}
},{"./selection/index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/create.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _creator = _interopRequireDefault(require("./creator"));

var _select = _interopRequireDefault(require("./select"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(name) {
  return (0, _select.default)((0, _creator.default)(name).call(document.documentElement));
}
},{"./creator":"../node_modules/d3-selection/src/creator.js","./select":"../node_modules/d3-selection/src/select.js"}],"../node_modules/d3-selection/src/local.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = local;
var nextId = 0;

function local() {
  return new Local();
}

function Local() {
  this._ = "@" + (++nextId).toString(36);
}

Local.prototype = local.prototype = {
  constructor: Local,
  get: function (node) {
    var id = this._;

    while (!(id in node)) if (!(node = node.parentNode)) return;

    return node[id];
  },
  set: function (node, value) {
    return node[this._] = value;
  },
  remove: function (node) {
    return this._ in node && delete node[this._];
  },
  toString: function () {
    return this._;
  }
};
},{}],"../node_modules/d3-selection/src/sourceEvent.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _on = require("./selection/on");

function _default() {
  var current = _on.event,
      source;

  while (source = current.sourceEvent) current = source;

  return current;
}
},{"./selection/on":"../node_modules/d3-selection/src/selection/on.js"}],"../node_modules/d3-selection/src/point.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

function _default(node, event) {
  var svg = node.ownerSVGElement || node;

  if (svg.createSVGPoint) {
    var point = svg.createSVGPoint();
    point.x = event.clientX, point.y = event.clientY;
    point = point.matrixTransform(node.getScreenCTM().inverse());
    return [point.x, point.y];
  }

  var rect = node.getBoundingClientRect();
  return [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
}
},{}],"../node_modules/d3-selection/src/mouse.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _sourceEvent = _interopRequireDefault(require("./sourceEvent"));

var _point = _interopRequireDefault(require("./point"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(node) {
  var event = (0, _sourceEvent.default)();
  if (event.changedTouches) event = event.changedTouches[0];
  return (0, _point.default)(node, event);
}
},{"./sourceEvent":"../node_modules/d3-selection/src/sourceEvent.js","./point":"../node_modules/d3-selection/src/point.js"}],"../node_modules/d3-selection/src/selectAll.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _index = require("./selection/index");

function _default(selector) {
  return typeof selector === "string" ? new _index.Selection([document.querySelectorAll(selector)], [document.documentElement]) : new _index.Selection([selector == null ? [] : selector], _index.root);
}
},{"./selection/index":"../node_modules/d3-selection/src/selection/index.js"}],"../node_modules/d3-selection/src/touch.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _sourceEvent = _interopRequireDefault(require("./sourceEvent"));

var _point = _interopRequireDefault(require("./point"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(node, touches, identifier) {
  if (arguments.length < 3) identifier = touches, touches = (0, _sourceEvent.default)().changedTouches;

  for (var i = 0, n = touches ? touches.length : 0, touch; i < n; ++i) {
    if ((touch = touches[i]).identifier === identifier) {
      return (0, _point.default)(node, touch);
    }
  }

  return null;
}
},{"./sourceEvent":"../node_modules/d3-selection/src/sourceEvent.js","./point":"../node_modules/d3-selection/src/point.js"}],"../node_modules/d3-selection/src/touches.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = _default;

var _sourceEvent = _interopRequireDefault(require("./sourceEvent"));

var _point = _interopRequireDefault(require("./point"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _default(node, touches) {
  if (touches == null) touches = (0, _sourceEvent.default)().touches;

  for (var i = 0, n = touches ? touches.length : 0, points = new Array(n); i < n; ++i) {
    points[i] = (0, _point.default)(node, touches[i]);
  }

  return points;
}
},{"./sourceEvent":"../node_modules/d3-selection/src/sourceEvent.js","./point":"../node_modules/d3-selection/src/point.js"}],"../node_modules/d3-selection/src/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
Object.defineProperty(exports, "create", {
  enumerable: true,
  get: function () {
    return _create.default;
  }
});
Object.defineProperty(exports, "creator", {
  enumerable: true,
  get: function () {
    return _creator.default;
  }
});
Object.defineProperty(exports, "local", {
  enumerable: true,
  get: function () {
    return _local.default;
  }
});
Object.defineProperty(exports, "matcher", {
  enumerable: true,
  get: function () {
    return _matcher.default;
  }
});
Object.defineProperty(exports, "mouse", {
  enumerable: true,
  get: function () {
    return _mouse.default;
  }
});
Object.defineProperty(exports, "namespace", {
  enumerable: true,
  get: function () {
    return _namespace.default;
  }
});
Object.defineProperty(exports, "namespaces", {
  enumerable: true,
  get: function () {
    return _namespaces.default;
  }
});
Object.defineProperty(exports, "clientPoint", {
  enumerable: true,
  get: function () {
    return _point.default;
  }
});
Object.defineProperty(exports, "select", {
  enumerable: true,
  get: function () {
    return _select.default;
  }
});
Object.defineProperty(exports, "selectAll", {
  enumerable: true,
  get: function () {
    return _selectAll.default;
  }
});
Object.defineProperty(exports, "selection", {
  enumerable: true,
  get: function () {
    return _index.default;
  }
});
Object.defineProperty(exports, "selector", {
  enumerable: true,
  get: function () {
    return _selector.default;
  }
});
Object.defineProperty(exports, "selectorAll", {
  enumerable: true,
  get: function () {
    return _selectorAll.default;
  }
});
Object.defineProperty(exports, "style", {
  enumerable: true,
  get: function () {
    return _style.styleValue;
  }
});
Object.defineProperty(exports, "touch", {
  enumerable: true,
  get: function () {
    return _touch.default;
  }
});
Object.defineProperty(exports, "touches", {
  enumerable: true,
  get: function () {
    return _touches.default;
  }
});
Object.defineProperty(exports, "window", {
  enumerable: true,
  get: function () {
    return _window.default;
  }
});
Object.defineProperty(exports, "event", {
  enumerable: true,
  get: function () {
    return _on.event;
  }
});
Object.defineProperty(exports, "customEvent", {
  enumerable: true,
  get: function () {
    return _on.customEvent;
  }
});

var _create = _interopRequireDefault(require("./create"));

var _creator = _interopRequireDefault(require("./creator"));

var _local = _interopRequireDefault(require("./local"));

var _matcher = _interopRequireDefault(require("./matcher"));

var _mouse = _interopRequireDefault(require("./mouse"));

var _namespace = _interopRequireDefault(require("./namespace"));

var _namespaces = _interopRequireDefault(require("./namespaces"));

var _point = _interopRequireDefault(require("./point"));

var _select = _interopRequireDefault(require("./select"));

var _selectAll = _interopRequireDefault(require("./selectAll"));

var _index = _interopRequireDefault(require("./selection/index"));

var _selector = _interopRequireDefault(require("./selector"));

var _selectorAll = _interopRequireDefault(require("./selectorAll"));

var _style = require("./selection/style");

var _touch = _interopRequireDefault(require("./touch"));

var _touches = _interopRequireDefault(require("./touches"));

var _window = _interopRequireDefault(require("./window"));

var _on = require("./selection/on");

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }
},{"./create":"../node_modules/d3-selection/src/create.js","./creator":"../node_modules/d3-selection/src/creator.js","./local":"../node_modules/d3-selection/src/local.js","./matcher":"../node_modules/d3-selection/src/matcher.js","./mouse":"../node_modules/d3-selection/src/mouse.js","./namespace":"../node_modules/d3-selection/src/namespace.js","./namespaces":"../node_modules/d3-selection/src/namespaces.js","./point":"../node_modules/d3-selection/src/point.js","./select":"../node_modules/d3-selection/src/select.js","./selectAll":"../node_modules/d3-selection/src/selectAll.js","./selection/index":"../node_modules/d3-selection/src/selection/index.js","./selector":"../node_modules/d3-selection/src/selector.js","./selectorAll":"../node_modules/d3-selection/src/selectorAll.js","./selection/style":"../node_modules/d3-selection/src/selection/style.js","./touch":"../node_modules/d3-selection/src/touch.js","./touches":"../node_modules/d3-selection/src/touches.js","./window":"../node_modules/d3-selection/src/window.js","./selection/on":"../node_modules/d3-selection/src/selection/on.js"}],"assets/world-110m.json":[function(require,module,exports) {
module.exports = {
  "type": "Topology",
  "objects": {
    "countries": {
      "type": "GeometryCollection",
      "geometries": [{
        "type": "Polygon",
        "arcs": [[0, 1, 2, 3, 4, 5]],
        "id": "004"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[6, 7, 8, 9]], [[10, 11, 12]]],
        "id": "024"
      }, {
        "type": "Polygon",
        "arcs": [[13, 14, 15, 16, 17]],
        "id": "008"
      }, {
        "type": "Polygon",
        "arcs": [[18, 19, 20, 21, 22]],
        "id": "784"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[23, 24]], [[25, 26, 27, 28, 29, 30]]],
        "id": "032"
      }, {
        "type": "Polygon",
        "arcs": [[31, 32, 33, 34, 35]],
        "id": "051"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[36]], [[37]], [[38]], [[39]], [[40]], [[41]], [[42]], [[43]]],
        "id": "010"
      }, {
        "type": "Polygon",
        "arcs": [[44]],
        "id": "260"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[45]], [[46]]],
        "id": "036"
      }, {
        "type": "Polygon",
        "arcs": [[47, 48, 49, 50, 51, 52, 53]],
        "id": "040"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[54, -35]], [[55, 56, -33, 57, 58]]],
        "id": "031"
      }, {
        "type": "Polygon",
        "arcs": [[59, 60, 61]],
        "id": "108"
      }, {
        "type": "Polygon",
        "arcs": [[62, 63, 64, 65, 66]],
        "id": "056"
      }, {
        "type": "Polygon",
        "arcs": [[67, 68, 69, 70, 71]],
        "id": "204"
      }, {
        "type": "Polygon",
        "arcs": [[72, 73, 74, -70, 75, 76]],
        "id": "854"
      }, {
        "type": "Polygon",
        "arcs": [[77, 78, 79]],
        "id": "050"
      }, {
        "type": "Polygon",
        "arcs": [[80, 81, 82, 83, 84, 85]],
        "id": "100"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[86]], [[87]], [[88]]],
        "id": "044"
      }, {
        "type": "Polygon",
        "arcs": [[89, 90, 91]],
        "id": "070"
      }, {
        "type": "Polygon",
        "arcs": [[92, 93, 94, 95, 96]],
        "id": "112"
      }, {
        "type": "Polygon",
        "arcs": [[97, 98, 99]],
        "id": "084"
      }, {
        "type": "Polygon",
        "arcs": [[100, 101, 102, 103, -31]],
        "id": "068"
      }, {
        "type": "Polygon",
        "arcs": [[-27, 104, -103, 105, 106, 107, 108, 109, 110, 111, 112]],
        "id": "076"
      }, {
        "type": "Polygon",
        "arcs": [[113, 114]],
        "id": "096"
      }, {
        "type": "Polygon",
        "arcs": [[115, 116]],
        "id": "064"
      }, {
        "type": "Polygon",
        "arcs": [[117, 118, 119, 120]],
        "id": "072"
      }, {
        "type": "Polygon",
        "arcs": [[121, 122, 123, 124, 125, 126, 127]],
        "id": "140"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[128]], [[129]], [[130]], [[131]], [[132]], [[133]], [[134]], [[135]], [[136]], [[137]], [[138, 139, 140, 141]], [[142]], [[143]], [[144]], [[145]], [[146]], [[147]], [[148]], [[149]], [[150]], [[151]], [[152]], [[153]], [[154]], [[155]], [[156]], [[157]], [[158]], [[159]], [[160]]],
        "id": "124"
      }, {
        "type": "Polygon",
        "arcs": [[-51, 161, 162, 163]],
        "id": "756"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[-24, 164]], [[-30, 165, 166, -101]]],
        "id": "152"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[167]], [[168, 169, 170, 171, 172, 173, -117, 174, 175, 176, 177, -4, 178, 179, 180, 181, 182, 183]]],
        "id": "156"
      }, {
        "type": "Polygon",
        "arcs": [[184, 185, 186, 187, -73, 188]],
        "id": "384"
      }, {
        "type": "Polygon",
        "arcs": [[189, 190, 191, 192, 193, 194, -128, 195]],
        "id": "120"
      }, {
        "type": "Polygon",
        "arcs": [[196, 197, -60, 198, 199, 200, 201, -10, 202, -13, 203, -126, 204]],
        "id": "180"
      }, {
        "type": "Polygon",
        "arcs": [[-12, 205, 206, -196, -127, -204]],
        "id": "178"
      }, {
        "type": "Polygon",
        "arcs": [[207, 208, 209, 210, 211, -107, 212]],
        "id": "170"
      }, {
        "type": "Polygon",
        "arcs": [[213, 214, 215, 216]],
        "id": "188"
      }, {
        "type": "Polygon",
        "arcs": [[217]],
        "id": "192"
      }, {
        "type": "Polygon",
        "arcs": [[218, 219]],
        "id": "-99"
      }, {
        "type": "Polygon",
        "arcs": [[220, -220]],
        "id": "196"
      }, {
        "type": "Polygon",
        "arcs": [[-53, 221, 222, 223]],
        "id": "203"
      }, {
        "type": "Polygon",
        "arcs": [[224, 225, -222, -52, -164, 226, 227, -64, 228, 229, 230]],
        "id": "276"
      }, {
        "type": "Polygon",
        "arcs": [[231, 232, 233, 234]],
        "id": "262"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[235]], [[-231, 236]]],
        "id": "208"
      }, {
        "type": "Polygon",
        "arcs": [[237, 238]],
        "id": "214"
      }, {
        "type": "Polygon",
        "arcs": [[239, 240, 241, 242, 243, 244, 245, 246]],
        "id": "012"
      }, {
        "type": "Polygon",
        "arcs": [[247, -208, 248]],
        "id": "218"
      }, {
        "type": "Polygon",
        "arcs": [[249, 250, 251, 252, 253]],
        "id": "818"
      }, {
        "type": "Polygon",
        "arcs": [[254, 255, 256, -235]],
        "id": "232"
      }, {
        "type": "Polygon",
        "arcs": [[257, 258, 259, 260]],
        "id": "724"
      }, {
        "type": "Polygon",
        "arcs": [[261, 262, 263]],
        "id": "233"
      }, {
        "type": "Polygon",
        "arcs": [[-234, 264, 265, 266, 267, 268, 269, -255]],
        "id": "231"
      }, {
        "type": "Polygon",
        "arcs": [[270, 271, 272, 273]],
        "id": "246"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[274]], [[275]]],
        "id": "242"
      }, {
        "type": "Polygon",
        "arcs": [[276]],
        "id": "238"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[277, 278, 279, -111]], [[280]], [[281, -227, -163, 282, 283, -259, 284, -66]]],
        "id": "250"
      }, {
        "type": "Polygon",
        "arcs": [[285, 286, -190, -207]],
        "id": "266"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[287, 288]], [[289]]],
        "id": "826"
      }, {
        "type": "Polygon",
        "arcs": [[290, 291, -58, -32, 292]],
        "id": "268"
      }, {
        "type": "Polygon",
        "arcs": [[293, -189, -77, 294]],
        "id": "288"
      }, {
        "type": "Polygon",
        "arcs": [[295, 296, 297, 298, 299, 300, -187]],
        "id": "324"
      }, {
        "type": "Polygon",
        "arcs": [[301, 302]],
        "id": "270"
      }, {
        "type": "Polygon",
        "arcs": [[303, 304, -299]],
        "id": "624"
      }, {
        "type": "Polygon",
        "arcs": [[305, -191, -287]],
        "id": "226"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[306]], [[307, -15, 308, -84, 309]]],
        "id": "300"
      }, {
        "type": "Polygon",
        "arcs": [[310]],
        "id": "304"
      }, {
        "type": "Polygon",
        "arcs": [[311, 312, -100, 313, 314, 315]],
        "id": "320"
      }, {
        "type": "Polygon",
        "arcs": [[316, 317, -109, 318]],
        "id": "328"
      }, {
        "type": "Polygon",
        "arcs": [[319, 320, -315, 321, 322]],
        "id": "340"
      }, {
        "type": "Polygon",
        "arcs": [[323, -92, 324, 325, 326, 327]],
        "id": "191"
      }, {
        "type": "Polygon",
        "arcs": [[-239, 328]],
        "id": "332"
      }, {
        "type": "Polygon",
        "arcs": [[-48, 329, 330, 331, 332, -328, 333]],
        "id": "348"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[334]], [[335, 336]], [[337]], [[338]], [[339]], [[340]], [[341]], [[342]], [[343, 344]], [[345]], [[346]], [[347, 348]], [[349]]],
        "id": "360"
      }, {
        "type": "Polygon",
        "arcs": [[-177, 350, -175, -116, -174, 351, -80, 352, 353]],
        "id": "356"
      }, {
        "type": "Polygon",
        "arcs": [[354, -288]],
        "id": "372"
      }, {
        "type": "Polygon",
        "arcs": [[355, -6, 356, 357, 358, 359, -55, -34, -57, 360]],
        "id": "364"
      }, {
        "type": "Polygon",
        "arcs": [[361, 362, 363, 364, 365, 366, -359]],
        "id": "368"
      }, {
        "type": "Polygon",
        "arcs": [[367]],
        "id": "352"
      }, {
        "type": "Polygon",
        "arcs": [[368, 369, 370, -254, 371, 372, 373]],
        "id": "376"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[374]], [[375]], [[376, 377, -283, -162, -50]]],
        "id": "380"
      }, {
        "type": "Polygon",
        "arcs": [[378]],
        "id": "388"
      }, {
        "type": "Polygon",
        "arcs": [[-369, 379, -365, 380, 381, -371, 382]],
        "id": "400"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[383]], [[384]], [[385]]],
        "id": "392"
      }, {
        "type": "Polygon",
        "arcs": [[386, 387, 388, 389, -181, 390]],
        "id": "398"
      }, {
        "type": "Polygon",
        "arcs": [[391, 392, 393, 394, -267, 395]],
        "id": "404"
      }, {
        "type": "Polygon",
        "arcs": [[-391, -180, 396, 397]],
        "id": "417"
      }, {
        "type": "Polygon",
        "arcs": [[398, 399, 400, 401]],
        "id": "116"
      }, {
        "type": "Polygon",
        "arcs": [[402, 403]],
        "id": "410"
      }, {
        "type": "Polygon",
        "arcs": [[-18, 404, 405, 406]],
        "id": "-99"
      }, {
        "type": "Polygon",
        "arcs": [[407, 408, -363]],
        "id": "414"
      }, {
        "type": "Polygon",
        "arcs": [[409, 410, -172, 411, -400]],
        "id": "418"
      }, {
        "type": "Polygon",
        "arcs": [[-373, 412, 413]],
        "id": "422"
      }, {
        "type": "Polygon",
        "arcs": [[414, 415, -296, -186]],
        "id": "430"
      }, {
        "type": "Polygon",
        "arcs": [[416, -247, 417, 418, -252, 419, 420]],
        "id": "434"
      }, {
        "type": "Polygon",
        "arcs": [[421]],
        "id": "144"
      }, {
        "type": "Polygon",
        "arcs": [[422]],
        "id": "426"
      }, {
        "type": "Polygon",
        "arcs": [[423, 424, 425, -93, 426]],
        "id": "440"
      }, {
        "type": "Polygon",
        "arcs": [[-228, -282, -65]],
        "id": "442"
      }, {
        "type": "Polygon",
        "arcs": [[427, -264, 428, -94, -426]],
        "id": "428"
      }, {
        "type": "Polygon",
        "arcs": [[-244, 429, 430]],
        "id": "504"
      }, {
        "type": "Polygon",
        "arcs": [[431, 432]],
        "id": "498"
      }, {
        "type": "Polygon",
        "arcs": [[433]],
        "id": "450"
      }, {
        "type": "Polygon",
        "arcs": [[434, -98, -313, 435, 436]],
        "id": "484"
      }, {
        "type": "Polygon",
        "arcs": [[-407, 437, -85, -309, -14]],
        "id": "807"
      }, {
        "type": "Polygon",
        "arcs": [[438, -241, 439, -74, -188, -301, 440]],
        "id": "466"
      }, {
        "type": "Polygon",
        "arcs": [[441, -78, -352, -173, -411, 442]],
        "id": "104"
      }, {
        "type": "Polygon",
        "arcs": [[443, -325, -91, 444, -405, -17]],
        "id": "499"
      }, {
        "type": "Polygon",
        "arcs": [[445, -183]],
        "id": "496"
      }, {
        "type": "Polygon",
        "arcs": [[446, 447, 448, 449, 450, 451, 452, 453]],
        "id": "508"
      }, {
        "type": "Polygon",
        "arcs": [[454, 455, 456, -242, -439]],
        "id": "478"
      }, {
        "type": "Polygon",
        "arcs": [[-454, 457, 458]],
        "id": "454"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[459, 460]], [[-348, 461, -115, 462]]],
        "id": "458"
      }, {
        "type": "Polygon",
        "arcs": [[463, -8, 464, -119, 465]],
        "id": "516"
      }, {
        "type": "Polygon",
        "arcs": [[466]],
        "id": "540"
      }, {
        "type": "Polygon",
        "arcs": [[-75, -440, -240, -417, 467, -194, 468, -71]],
        "id": "562"
      }, {
        "type": "Polygon",
        "arcs": [[469, -72, -469, -193]],
        "id": "566"
      }, {
        "type": "Polygon",
        "arcs": [[470, -323, 471, -215]],
        "id": "558"
      }, {
        "type": "Polygon",
        "arcs": [[-229, -63, 472]],
        "id": "528"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[473, -274, 474, 475]], [[476]], [[477]], [[478]]],
        "id": "578"
      }, {
        "type": "Polygon",
        "arcs": [[-351, -176]],
        "id": "524"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[479]], [[480]]],
        "id": "554"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[481, 482, -22, 483]], [[-20, 484]]],
        "id": "512"
      }, {
        "type": "Polygon",
        "arcs": [[-178, -354, 485, -357, -5]],
        "id": "586"
      }, {
        "type": "Polygon",
        "arcs": [[486, -217, 487, -210]],
        "id": "591"
      }, {
        "type": "Polygon",
        "arcs": [[-167, 488, -249, -213, -106, -102]],
        "id": "604"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[489]], [[490]], [[491]], [[492]], [[493]], [[494]], [[495]]],
        "id": "608"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[496]], [[497]], [[-344, 498]], [[499]]],
        "id": "598"
      }, {
        "type": "Polygon",
        "arcs": [[-226, 500, 501, -427, -97, 502, 503, -223]],
        "id": "616"
      }, {
        "type": "Polygon",
        "arcs": [[504]],
        "id": "630"
      }, {
        "type": "Polygon",
        "arcs": [[505, 506, -404, 507, -169]],
        "id": "408"
      }, {
        "type": "Polygon",
        "arcs": [[-261, 508]],
        "id": "620"
      }, {
        "type": "Polygon",
        "arcs": [[-104, -105, -26]],
        "id": "600"
      }, {
        "type": "Polygon",
        "arcs": [[-383, -370]],
        "id": "275"
      }, {
        "type": "Polygon",
        "arcs": [[509, 510]],
        "id": "634"
      }, {
        "type": "Polygon",
        "arcs": [[511, -433, 512, 513, -81, 514, -332]],
        "id": "642"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[515]], [[-502, 516, -424]], [[517]], [[518]], [[519]], [[520]], [[521]], [[-506, -184, -446, -182, -390, 522, -59, -292, 523, 524, -95, -429, -263, 525, -271, -474, 526]], [[527]], [[528]], [[529]]],
        "id": "643"
      }, {
        "type": "Polygon",
        "arcs": [[530, -61, -198, 531]],
        "id": "646"
      }, {
        "type": "Polygon",
        "arcs": [[-243, -457, 532, -430]],
        "id": "732"
      }, {
        "type": "Polygon",
        "arcs": [[533, -381, -364, -409, 534, -511, 535, -23, -483, 536]],
        "id": "682"
      }, {
        "type": "Polygon",
        "arcs": [[537, 538, -123, 539, -420, -251, 540, -256, -270, 541]],
        "id": "729"
      }, {
        "type": "Polygon",
        "arcs": [[542, -268, -395, 543, -205, -125, 544, -538]],
        "id": "728"
      }, {
        "type": "Polygon",
        "arcs": [[545, -455, -441, -300, -305, 546, -303]],
        "id": "686"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[547]], [[548]], [[549]], [[550]], [[551]]],
        "id": "090"
      }, {
        "type": "Polygon",
        "arcs": [[552, -297, -416]],
        "id": "694"
      }, {
        "type": "Polygon",
        "arcs": [[553, -316, -321]],
        "id": "222"
      }, {
        "type": "Polygon",
        "arcs": [[-265, -233, 554, 555]],
        "id": "-99"
      }, {
        "type": "Polygon",
        "arcs": [[-396, -266, -556, 556]],
        "id": "706"
      }, {
        "type": "Polygon",
        "arcs": [[-86, -438, -406, -445, -90, -324, -333, -515]],
        "id": "688"
      }, {
        "type": "Polygon",
        "arcs": [[557, -279, 558, -110, -318]],
        "id": "740"
      }, {
        "type": "Polygon",
        "arcs": [[-504, 559, -330, -54, -224]],
        "id": "703"
      }, {
        "type": "Polygon",
        "arcs": [[-49, -334, -327, 560, -377]],
        "id": "705"
      }, {
        "type": "Polygon",
        "arcs": [[-475, -273, 561]],
        "id": "752"
      }, {
        "type": "Polygon",
        "arcs": [[562, -450]],
        "id": "748"
      }, {
        "type": "Polygon",
        "arcs": [[-380, -374, -414, 563, 564, -366]],
        "id": "760"
      }, {
        "type": "Polygon",
        "arcs": [[-468, -421, -540, -122, -195]],
        "id": "148"
      }, {
        "type": "Polygon",
        "arcs": [[565, -295, -76, -69]],
        "id": "768"
      }, {
        "type": "Polygon",
        "arcs": [[566, -461, 567, -443, -410, -399]],
        "id": "764"
      }, {
        "type": "Polygon",
        "arcs": [[-397, -179, -3, 568]],
        "id": "762"
      }, {
        "type": "Polygon",
        "arcs": [[-356, 569, -388, 570, -1]],
        "id": "795"
      }, {
        "type": "Polygon",
        "arcs": [[571, -336]],
        "id": "626"
      }, {
        "type": "Polygon",
        "arcs": [[572]],
        "id": "780"
      }, {
        "type": "Polygon",
        "arcs": [[-246, 573, -418]],
        "id": "788"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[-293, -36, -360, -367, -565, 574]], [[-310, -83, 575]]],
        "id": "792"
      }, {
        "type": "Polygon",
        "arcs": [[576]],
        "id": "158"
      }, {
        "type": "Polygon",
        "arcs": [[-393, 577, -447, -459, 578, -201, 579, -199, -62, -531, 580]],
        "id": "834"
      }, {
        "type": "Polygon",
        "arcs": [[-532, -197, -544, -394, -581]],
        "id": "800"
      }, {
        "type": "Polygon",
        "arcs": [[-525, 581, -513, -432, -512, -331, -560, -503, -96]],
        "id": "804"
      }, {
        "type": "Polygon",
        "arcs": [[-113, 582, -28]],
        "id": "858"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[583]], [[584]], [[585]], [[586]], [[587]], [[588, -437, 589, -139]], [[590]], [[591]], [[592]], [[-141, 593]]],
        "id": "840"
      }, {
        "type": "Polygon",
        "arcs": [[-571, -387, -398, -569, -2]],
        "id": "860"
      }, {
        "type": "Polygon",
        "arcs": [[594, -319, -108, -212]],
        "id": "862"
      }, {
        "type": "Polygon",
        "arcs": [[595, -401, -412, -171]],
        "id": "704"
      }, {
        "type": "MultiPolygon",
        "arcs": [[[596]], [[597]]],
        "id": "548"
      }, {
        "type": "Polygon",
        "arcs": [[598, -537, -482]],
        "id": "887"
      }, {
        "type": "Polygon",
        "arcs": [[-466, -118, 599, -451, -563, -449, 600], [-423]],
        "id": "710"
      }, {
        "type": "Polygon",
        "arcs": [[-458, -453, 601, -120, -465, -7, -202, -579]],
        "id": "894"
      }, {
        "type": "Polygon",
        "arcs": [[-600, -121, -602, -452]],
        "id": "716"
      }]
    },
    "land": {
      "type": "GeometryCollection",
      "geometries": [{
        "type": "MultiPolygon",
        "arcs": [[[595, 401, 566, 459, 567, 441, 78, 352, 485, 357, 361, 407, 534, 509, 535, 18, 484, 20, 483, 598, 533, 381, 249, 540, 256, 231, 554, 556, 391, 577, 447, 600, 463, 8, 202, 10, 205, 285, 305, 191, 469, 67, 565, 293, 184, 414, 552, 297, 303, 546, 301, 545, 455, 532, 430, 244, 573, 418, 252, 371, 412, 563, 574, 290, 523, 581, 513, 81, 575, 307, 15, 443, 325, 560, 377, 283, 259, 508, 257, 284, 66, 472, 229, 236, 224, 500, 516, 424, 427, 261, 525, 271, 561, 475, 526, 506, 402, 507, 169], [123, 544, 538], [199, 579], [542, 268, 541], [388, 522, 55, 360, 569]], [[24, 164]], [[582, 28, 165, 488, 247, 208, 486, 213, 470, 319, 553, 311, 435, 589, 139, 593, 141, 588, 434, 98, 313, 321, 471, 215, 487, 210, 594, 316, 557, 279, 111], [558, 277]], [[36]], [[37]], [[38]], [[39]], [[40]], [[41]], [[42]], [[43]], [[44]], [[45]], [[46]], [[86]], [[87]], [[88]], [[461, 113, 462, 348]], [[128]], [[129]], [[130]], [[131]], [[132]], [[133]], [[134]], [[135]], [[136]], [[137]], [[142]], [[143]], [[144]], [[145]], [[146]], [[147]], [[148]], [[149]], [[150]], [[151]], [[152]], [[153]], [[154]], [[155]], [[156]], [[157]], [[158]], [[159]], [[160]], [[167]], [[217]], [[218, 220]], [[235]], [[237, 328]], [[274]], [[275]], [[276]], [[280]], [[288, 354]], [[289]], [[306]], [[310]], [[334]], [[336, 571]], [[337]], [[338]], [[339]], [[340]], [[341]], [[342]], [[344, 498]], [[345]], [[346]], [[349]], [[367]], [[374]], [[375]], [[378]], [[383]], [[384]], [[385]], [[421]], [[433]], [[466]], [[476]], [[477]], [[478]], [[479]], [[480]], [[489]], [[490]], [[491]], [[492]], [[493]], [[494]], [[495]], [[496]], [[497]], [[499]], [[504]], [[515]], [[517]], [[518]], [[519]], [[520]], [[521]], [[527]], [[528]], [[529]], [[547]], [[548]], [[549]], [[550]], [[551]], [[572]], [[576]], [[583]], [[584]], [[585]], [[586]], [[587]], [[590]], [[591]], [[592]], [[596]], [[597]]]
      }]
    }
  },
  "arcs": [[[67002, 71642], [284, -224], [209, 79], [58, 268], [219, 89], [157, 180], [55, 472], [234, 114], [44, 211], [131, -158], [84, -19]], [[68477, 72654], [154, -4], [210, -124]], [[68841, 72526], [85, -72], [201, 189], [93, -114], [90, 271], [166, -12], [43, 86], [29, 239], [120, 205], [150, -134], [-30, -181], [84, -28], [-26, -496], [110, -194], [97, 125], [123, 58], [173, 265], [192, -44], [286, -1]], [[70827, 72688], [50, -169]], [[70877, 72519], [-162, -67], [-141, -109], [-319, -68], [-298, -124], [-163, -258], [66, -250], [32, -294], [-139, -248], [12, -227], [-76, -213], [-265, 18], [110, -390], [-177, -150], [-118, -356], [15, -355], [-108, -166], [-103, 55], [-212, -77], [-31, -166], [-207, 1], [-154, -334], [-10, -503], [-361, -246], [-194, 52], [-56, -129], [-166, 75], [-278, -88], [-465, 301]], [[66909, 68203], [252, 536], [-23, 380], [-210, 100], [-22, 375], [-91, 472], [119, 323], [-121, 87], [76, 430], [113, 736]], [[56642, 44124], [29, -184], [-32, -286], [49, -277], [-41, -222], [24, -203], [-579, 7], [-13, -1880], [188, -483], [181, -369]], [[56448, 40227], [-510, -241], [-673, 83], [-192, 284], [-1126, -26], [-42, -41], [-166, 267], [-180, 17], [-166, -100], [-134, -113]], [[53259, 40357], [-26, 372], [38, 519], [96, 541], [15, 254], [90, 532], [66, 243], [159, 386], [90, 263], [29, 438], [-15, 335], [-83, 211], [-74, 358], [-68, 355], [15, 122], [85, 235], [-84, 570], [-57, 396], [-139, 374], [26, 115]], [[53422, 46976], [115, 79], [80, -11], [98, 71], [820, -8], [68, -440], [80, -354], [64, -191], [106, -309], [184, 47], [91, 83], [154, -83], [42, 148], [69, 344], [172, 23], [15, 103], [142, 2], [-24, -213], [337, 5], [5, -372], [56, -228], [-41, -356], [21, -363], [93, -219], [-15, -703], [68, 54], [121, -15], [172, 89], [127, -35]], [[53383, 47159], [-74, 444]], [[53309, 47603], [112, 255], [84, 100], [104, -203]], [[53609, 47755], [-101, -124], [-45, -152], [-9, -258], [-71, -62]], [[55719, 75309], [-35, -201], [39, -254], [115, -144]], [[55838, 74710], [-5, -155], [-91, -85], [-16, -192], [-129, -287]], [[55597, 73991], [-48, 41], [-5, 130], [-154, 199], [-24, 281], [23, 403], [38, 184], [-47, 93]], [[55380, 75322], [-18, 188], [120, 291], [18, -111], [75, 52]], [[55575, 75742], [59, -159], [66, -60], [19, -214]], [[64327, 64904], [49, 29], [11, -162], [217, 93], [230, -15], [168, -18], [190, 400], [207, 379], [176, 364]], [[65575, 65974], [52, -202]], [[65627, 65772], [38, -466]], [[65665, 65306], [-142, -3], [-23, -384], [50, -82], [-126, -117], [-1, -241], [-81, -245], [-7, -238]], [[65335, 63996], [-56, -125], [-835, 298], [-106, 599], [-11, 136]], [[31400, 18145], [-168, 16], [-297, 1], [0, 1319]], [[30935, 19481], [106, -274], [139, -443], [361, -355], [389, -147], [-125, -296], [-264, -29], [-141, 208]], [[32587, 37434], [511, -964], [227, -89], [339, -437], [286, -231], [40, -261], [-273, -898], [280, -160], [312, -91], [220, 95], [252, 453], [45, 521]], [[34826, 35372], [138, 114], [139, -341], [-6, -472], [-234, -326], [-186, -241], [-314, -573], [-370, -806]], [[33993, 32727], [-70, -473], [-74, -607], [3, -588], [-61, -132], [-21, -382]], [[33770, 30545], [-19, -308], [353, -506], [-38, -408], [173, -257], [-14, -289], [-267, -757], [-412, -317], [-557, -123], [-305, 59], [59, -352], [-57, -442], [51, -298], [-167, -208], [-284, -82], [-267, 216], [-108, -155], [39, -587], [188, -178], [152, 186], [82, -307], [-255, -183], [-223, -367], [-41, -595], [-66, -316], [-262, -2], [-218, -302], [-80, -443], [273, -433], [266, -119], [-96, -531], [-328, -333], [-180, -692], [-254, -234], [-113, -276], [89, -614], [185, -342], [-117, 30]], [[30952, 19680], [-257, 93], [-672, 79], [-115, 344], [6, 443], [-185, -38], [-98, 214], [-24, 626], [213, 260], [88, 375], [-33, 299], [148, 504], [101, 782], [-30, 347], [122, 112], [-30, 223], [-129, 118], [92, 248], [-126, 224], [-65, 682], [112, 120], [-47, 720], [65, 605], [75, 527], [166, 215], [-84, 576], [-1, 543], [210, 386], [-7, 494], [159, 576], [1, 544], [-72, 108], [-128, 1020], [171, 607], [-27, 572], [100, 537], [182, 555], [196, 367], [-83, 232], [58, 190], [-9, 985], [302, 291], [96, 614], [-34, 148]], [[31359, 37147], [231, 534], [364, -144], [163, -427], [109, 475], [316, -24], [45, -127]], [[62106, 74858], [386, 92]], [[62492, 74950], [57, -155], [106, -103], [-56, -148], [148, -202], [-78, -189], [118, -160], [124, -97], [7, -410]], [[62918, 73486], [-101, -17]], [[62817, 73469], [-113, 342], [1, 91], [-123, -2], [-82, 159], [-58, -16]], [[62442, 74043], [-109, 172], [-207, 147], [27, 288], [-47, 208]], [[33452, 3290], [-82, -301], [-81, -266], [-582, 81], [-621, -35], [-348, 197], [0, 23], [-152, 174], [625, -23], [599, -58], [207, 243], [147, 208], [288, -243]], [[5775, 3611], [-533, -81], [-364, 208], [-163, 209], [-11, 35], [-180, 162], [169, 220], [517, -93], [277, -185], [212, -209], [76, -266]], [[37457, 4468], [342, -255], [120, -359], [33, -254], [11, -301], [-430, -186], [-452, -150], [-522, -139], [-582, -116], [-658, 35], [-365, 197], [49, 243], [593, 162], [239, 197], [174, 254], [126, 220], [168, 209], [180, 243], [141, 0], [414, 127], [419, -127]], [[16330, 7154], [359, -93], [332, 104], [-158, -208], [-261, -151], [-386, 47], [-278, 208], [60, 197], [332, -104]], [[15122, 7165], [425, -231], [-164, 23], [-359, 58], [-381, 162], [202, 127], [277, -139]], [[22505, 8080], [305, -81], [304, 69], [163, -335], [-217, 46], [-337, -23], [-343, 23], [-376, -35], [-283, 116], [-146, 243], [174, 104], [353, -81], [403, -46]], [[30985, 8657], [33, -266], [-49, -231], [-76, -220], [-326, -81], [-311, -116], [-364, 11], [136, 232], [-327, -81], [-310, -81], [-212, 174], [-16, 243], [305, 231], [190, 70], [321, -23], [82, 301], [16, 219], [-6, 475], [158, 278], [256, 93], [147, -220], [65, -220], [120, -267], [92, -254], [76, -267]], [[0, 529], [16, -5], [245, 344], [501, -185], [32, 21], [294, 188], [38, -7], [32, -4], [402, -246], [352, 246], [63, 34], [816, 104], [265, -138], [130, -71], [419, -196], [789, -151], [625, -185], [1072, -139], [800, 162], [1181, -116], [669, -185], [734, 174], [773, 162], [60, 278], [-1094, 23], [-898, 139], [-234, 231], [-745, 128], [49, 266], [103, 243], [104, 220], [-55, 243], [-462, 162], [-212, 209], [-430, 185], [675, -35], [642, 93], [402, -197], [495, 173], [457, 220], [223, 197], [-98, 243], [-359, 162], [-408, 174], [-571, 35], [-500, 81], [-539, 58], [-180, 220], [-359, 185], [-217, 208], [-87, 672], [136, -58], [250, -185], [457, 58], [441, 81], [228, -255], [441, 58], [370, 127], [348, 162], [315, 197], [419, 58], [-11, 220], [-97, 220], [81, 208], [359, 104], [163, -196], [425, 115], [321, 151], [397, 12], [375, 57], [376, 139], [299, 128], [337, 127], [218, -35], [190, -46], [414, 81], [370, -104], [381, 11], [364, 81], [375, -57], [414, -58], [386, 23], [403, -12], [413, -11], [381, 23], [283, 174], [337, 92], [349, -127], [331, 104], [300, 208], [179, -185], [98, -208], [180, -197], [288, 174], [332, -220], [375, -70], [321, -162], [392, 35], [354, 104], [418, -23], [376, -81], [381, -104], [147, 254], [-180, 197], [-136, 209], [-359, 46], [-158, 220], [-60, 220], [-98, 440], [213, -81], [364, -35], [359, 35], [327, -93], [283, -174], [119, -208], [376, -35], [359, 81], [381, 116], [342, 70], [283, -139], [370, 46], [239, 451], [224, -266], [321, -104], [348, 58], [228, -232], [365, -23], [337, -69], [332, -128], [218, 220], [108, 209], [278, -232], [381, 58], [283, -127], [190, -197], [370, 58], [288, 127], [283, 151], [337, 81], [392, 69], [354, 81], [272, 127], [163, 186], [65, 254], [-32, 244], [-87, 231], [-98, 232], [-87, 231], [-71, 209], [-16, 231], [27, 232], [130, 220], [109, 243], [44, 231], [-55, 255], [-32, 232], [136, 266], [152, 173], [180, 220], [190, 186], [223, 173], [109, 255], [152, 162], [174, 151], [267, 34], [174, 186], [196, 115], [228, 70], [202, 150], [157, 186], [218, 69], [163, -151], [-103, -196], [-283, -174], [-120, -127], [-206, 92], [-229, -58], [-190, -139], [-202, -150], [-136, -174], [-38, -231], [17, -220], [130, -197], [-190, -139], [-261, -46], [-153, -197], [-163, -185], [-174, -255], [-44, -220], [98, -243], [147, -185], [229, -139], [212, -185], [114, -232], [60, -220], [82, -232], [130, -196], [82, -220], [38, -544], [81, -220], [22, -232], [87, -231], [-38, -313], [-152, -243], [-163, -197], [-370, -81], [-125, -208], [-169, -197], [-419, -220], [-370, -93], [-348, -127], [-376, -128], [-223, -243], [-446, -23], [-489, 23], [-441, -46], [-468, 0], [87, -232], [424, -104], [311, -162], [174, -208], [-310, -185], [-479, 58], [-397, -151], [-17, -243], [-11, -232], [327, -196], [60, -220], [353, -220], [588, -93], [500, -162], [398, -185], [506, -186], [690, -92], [681, -162], [473, -174], [517, -197], [272, -278], [136, -220], [337, 209], [457, 173], [484, 186], [577, 150], [495, 162], [691, 12], [680, -81], [560, -139], [180, 255], [386, 173], [702, 12], [550, 127], [522, 128], [577, 81], [614, 104], [430, 150], [-196, 209], [-119, 208], [0, 220], [-539, -23], [-571, -93], [-544, 0], [-77, 220], [39, 440], [125, 128], [397, 138], [468, 139], [337, 174], [337, 174], [251, 231], [380, 104], [376, 81], [190, 47], [430, 23], [408, 81], [343, 116], [337, 139], [305, 139], [386, 185], [245, 197], [261, 173], [82, 232], [-294, 139], [98, 243], [185, 185], [288, 116], [305, 139], [283, 185], [217, 232], [136, 277], [202, 163], [331, -35], [136, -197], [332, -23], [11, 220], [142, 231], [299, -58], [71, -220], [331, -34], [360, 104], [348, 69], [315, -34], [120, -243], [305, 196], [283, 105], [315, 81], [310, 81], [283, 139], [310, 92], [240, 128], [168, 208], [207, -151], [288, 81], [202, -277], [157, -209], [316, 116], [125, 232], [283, 162], [365, -35], [108, -220], [229, 220], [299, 69], [326, 23], [294, -11], [310, -70], [300, -34], [130, -197], [180, -174], [304, 104], [327, 24], [315, 0], [310, 11], [278, 81], [294, 70], [245, 162], [261, 104], [283, 58], [212, 162], [152, 324], [158, 197], [288, -93], [109, -208], [239, -139], [289, 46], [196, -208], [206, -151], [283, 139], [98, 255], [250, 104], [289, 197], [272, 81], [326, 116], [218, 127], [228, 139], [218, 127], [261, -69], [250, 208], [180, 162], [261, -11], [229, 139], [54, 208], [234, 162], [228, 116], [278, 93], [256, 46], [244, -35], [262, -58], [223, -162], [27, -254], [245, -197], [168, -162], [332, -70], [185, -162], [229, -162], [266, -35], [223, 116], [240, 243], [261, -127], [272, -70], [261, -69], [272, -46], [277, 0], [229, -614], [-11, -150], [-33, -267], [-266, -150], [-218, -220], [38, -232], [310, 12], [-38, -232], [-141, -220], [-131, -243], [212, -185], [321, -58], [321, 104], [153, 232], [92, 220], [153, 185], [174, 174], [70, 208], [147, 289], [174, 58], [316, 24], [277, 69], [283, 93], [136, 231], [82, 220], [190, 220], [272, 151], [234, 115], [153, 197], [157, 104], [202, 93], [277, -58], [250, 58], [272, 69], [305, -34], [201, 162], [142, 393], [103, -162], [131, -278], [234, -115], [266, -47], [267, 70], [283, -46], [261, -12], [174, 58], [234, -35], [212, -127], [250, 81], [300, 0], [255, 81], [289, -81], [185, 197], [141, 196], [191, 163], [348, 439], [179, -81], [212, -162], [185, -208], [354, -359], [272, -12], [256, 0], [299, 70], [299, 81], [229, 162], [190, 174], [310, 23], [207, 127], [218, -116], [141, -185], [196, -185], [305, 23], [190, -150], [332, -151], [348, -58], [288, 47], [218, 185], [185, 185], [250, 46], [251, -81], [288, -58], [261, 93], [250, 0], [245, -58], [256, -58], [250, 104], [299, 93], [283, 23], [316, 0], [255, 58], [251, 46], [76, 290], [11, 243], [174, -162], [49, -266], [92, -244], [115, -196], [234, -105], [315, 35], [365, 12], [250, 35], [364, 0], [262, 11], [364, -23], [310, -46], [196, -186], [-54, -220], [179, -173], [299, -139], [310, -151], [360, -104], [375, -92], [283, -93], [315, -12], [180, 197], [245, -162], [212, -185], [245, -139], [337, -58], [321, -69], [136, -232], [316, -139], [212, -208], [310, -93], [321, 12], [299, -35], [332, 12], [332, -47], [310, -81], [288, -139], [289, -116], [195, -173], [-32, -232], [-147, -208], [-125, -266], [-98, -209], [-131, -243], [-364, -93], [-163, -208], [-360, -127], [-125, -232], [-190, -220], [-201, -185], [-115, -243], [-70, -220], [-28, -266], [6, -220], [158, -232], [60, -220], [130, -208], [517, -81], [109, -255], [-501, -93], [-424, -127], [-528, -23], [-234, -336], [-49, -278], [-119, -220], [-147, -220], [370, -196], [141, -244], [239, -219], [338, -197], [386, -186], [419, -185], [636, -185], [142, -289], [800, -128], [53, -45], [208, -175], [767, 151], [636, -186], [-99520, -142]], [[69148, 21851], [179, -186], [263, -74], [9, -112], [-77, -269], [-427, -38], [-7, 314], [41, 244], [19, 121]], [[90387, 26479], [269, -204], [151, 81], [217, 113], [166, -39], [20, -702], [-95, -203], [-29, -476], [-97, 162], [-193, -412], [-57, 32], [-171, 19], [-171, 505], [-38, 390], [-160, 515], [7, 271], [181, -52]], [[89877, 42448], [100, -464], [179, 223], [92, -250], [133, -231], [-29, -262], [60, -506], [42, -295], [70, -72], [75, -505], [-27, -307], [90, -400], [301, -309], [197, -281], [186, -257], [-37, -143], [159, -371], [108, -639], [111, 130], [113, -256], [68, 91], [48, -626], [197, -363], [129, -226], [217, -478], [78, -475], [7, -337], [-19, -365], [132, -502], [-16, -523], [-48, -274], [-75, -527], [6, -339], [-55, -423], [-123, -538], [-205, -290], [-102, -458], [-93, -292], [-82, -510], [-107, -294], [-70, -442], [-36, -407], [14, -187], [-159, -205], [-311, -22], [-257, -242], [-127, -229], [-168, -254], [-230, 262], [-170, 104], [43, 308], [-152, -112], [-243, -428], [-240, 160], [-158, 94], [-159, 42], [-269, 171], [-179, 364], [-52, 449], [-64, 298], [-137, 240], [-267, 71], [91, 287], [-67, 438], [-136, -408], [-247, -109], [146, 327], [42, 341], [107, 289], [-22, 438], [-226, -504], [-174, -202], [-106, -470], [-217, 243], [9, 313], [-174, 429], [-147, 221], [52, 137], [-356, 358], [-195, 17], [-267, 287], [-498, -56], [-359, -211], [-317, -197], [-265, 39], [-294, -303], [-241, -137], [-53, -309], [-103, -240], [-236, -15], [-174, -52], [-246, 107], [-199, -64], [-191, -27], [-165, -315], [-81, 26], [-140, -167], [-133, -187], [-203, 23], [-186, 0], [-295, 377], [-149, 113], [6, 338], [138, 81], [47, 134], [-10, 212], [34, 411], [-31, 350], [-147, 598], [-45, 337], [12, 336], [-111, 385], [-7, 174], [-123, 235], [-35, 463], [-158, 467], [-39, 252], [122, -255], [-93, 548], [137, -171], [83, -229], [-5, 303], [-138, 465], [-26, 186], [-65, 177], [31, 341], [56, 146], [38, 295], [-29, 346], [114, 425], [21, -450], [118, 406], [225, 198], [136, 252], [212, 217], [126, 46], [77, -73], [219, 220], [168, 66], [42, 129], [74, 54], [153, -14], [292, 173], [151, 262], [71, 316], [163, 300], [13, 236], [7, 321], [194, 502], [117, -510], [119, 118], [-99, 279], [87, 287], [122, -128], [34, 449], [152, 291], [67, 233], [140, 101], [4, 165], [122, -69], [5, 148], [122, 85], [134, 80], [205, -271], [155, -350], [173, -4], [177, -56], [-59, 325], [133, 473], [126, 155], [-44, 147], [121, 338], [168, 208], [142, -70], [234, 111], [-5, 302], [-204, 195], [148, 86], [184, -147], [148, -242], [234, -151], [79, 60], [172, -182], [162, 169], [105, -51], [65, 113], [127, -292], [-74, -316], [-105, -239], [-96, -20], [32, -236], [-81, -295], [-99, -291], [20, -166], [221, -327], [214, -189], [143, -204], [201, -350], [78, 1], [145, -151], [43, -183], [265, -200], [183, 202], [55, 317], [56, 262], [34, 324], [85, 470], [-39, 286], [20, 171], [-32, 339], [37, 445], [53, 120], [-43, 197], [67, 313], [52, 325], [7, 168], [104, 222], [78, -289], [19, -371], [70, -71], [11, -249], [101, -300], [21, -335], [-10, -214]], [[54716, 79012], [-21, -241], [-156, -2], [53, -128], [-92, -380]], [[54500, 78261], [-53, -100], [-243, -14], [-140, -134], [-229, 45]], [[53835, 78058], [-398, 153], [-62, 205], [-274, -102], [-32, -113], [-169, 84]], [[52900, 78285], [-142, 16], [-125, 108], [42, 145], [-10, 104]], [[52665, 78658], [83, 33], [141, -164], [39, 156], [245, -25], [199, 106], [133, -18], [87, -121], [26, 100], [-40, 385], [100, 75], [98, 272]], [[53776, 79457], [206, -190], [157, 242], [98, 44], [215, -180], [131, 30], [128, -111]], [[54711, 79292], [-23, -75], [28, -205]], [[62817, 73469], [-190, 78], [-141, 273], [-44, 223]], [[63495, 75281], [146, -311], [141, -419], [130, -28], [85, -159], [-228, -47], [-49, -459], [-48, -207], [-101, -138], [7, -293]], [[63578, 73220], [-69, -29], [-173, 309], [95, 292], [-82, 174], [-104, -44], [-327, -436]], [[62492, 74950], [68, 96], [207, -169], [149, -36], [38, 70], [-136, 319], [72, 82]], [[62890, 75312], [78, -20], [191, -359], [122, -40], [48, 150], [166, 238]], [[58149, 47921], [-17, 713], [-70, 268]], [[58062, 48902], [169, -46], [85, 336], [147, -38]], [[58463, 49154], [16, -233], [60, -134], [3, -192], [-69, -124], [-108, -308], [-101, -214], [-115, -28]], [[50920, 80916], [204, -47], [257, 123], [176, -258], [153, -138]], [[51710, 80596], [-32, -400]], [[51678, 80196], [-72, -22], [-30, -331]], [[51576, 79843], [-243, 269], [-143, -46], [-194, 279], [-129, 237], [-129, 10], [-40, 207]], [[50698, 80799], [222, 117]], [[50747, 54278], [-229, -69]], [[50518, 54209], [-69, 407], [13, 1357], [-56, 122], [-11, 290], [-96, 207], [-85, 174], [35, 311]], [[50249, 57077], [96, 67], [56, 258], [136, 56], [61, 176]], [[50598, 57634], [93, 173], [100, 2], [212, -340]], [[51003, 57469], [-11, -197], [62, -350], [-54, -238], [29, -159], [-135, -366], [-86, -181], [-52, -372], [7, -376], [-16, -952]], [[49214, 56277], [-190, 152], [-130, -22], [-97, -149], [-125, 125], [-49, 195], [-125, 129]], [[48498, 56707], [-18, 343], [76, 250], [-7, 200], [221, 490], [41, 405], [76, 144], [134, -79], [116, 120], [38, 152], [216, 265], [53, 184], [259, 246], [153, 84], [70, -114], [178, 3]], [[50104, 59400], [-22, -286], [37, -269], [156, -386], [9, -286], [320, -134], [-6, -405]], [[50249, 57077], [-243, 13]], [[50006, 57090], [-128, 47], [-90, -96], [-123, 43], [-482, -27], [-7, -336], [38, -444]], [[75742, 63602], [-6, -424], [-97, 90], [18, -476]], [[75657, 62792], [-79, 308], [-16, 301], [-53, 285], [-116, 344], [-256, 23], [25, -243], [-87, -329], [-118, 120], [-41, -108], [-78, 65], [-108, 53]], [[74730, 63611], [-43, 486], [-96, 444], [47, 356], [-171, 159], [62, 215], [173, 220], [-200, 313], [98, 401], [220, -255], [133, -30], [24, -410], [265, -81], [257, 8], [160, -101], [-128, -500], [-124, -34], [-86, -336], [152, -306], [46, 377], [76, 2], [147, -937]], [[56293, 76715], [80, -243], [108, 43], [213, -92], [408, -31], [138, 150], [327, 138], [202, -215], [163, -62]], [[57932, 76403], [-144, -245], [-101, -422], [89, -337]], [[57776, 75399], [-239, 79], [-283, -186]], [[57254, 75292], [-3, -294], [-252, -56], [-196, 206], [-222, -162], [-206, 17]], [[56375, 75003], [-20, 391], [-139, 189]], [[56216, 75583], [46, 84], [-30, 70], [47, 188], [105, 185], [-135, 255], [-24, 216], [68, 134]], [[28462, 64617], [-68, -29], [-70, 340], [-104, 171], [60, 375], [84, -23], [97, -491], [1, -343]], [[28383, 66284], [-303, -95], [-19, 219], [130, 47], [184, -18], [8, -153]], [[28611, 66290], [-48, -420], [-51, 75], [4, 309], [-124, 234], [-1, 67], [220, -265]], [[55279, 77084], [100, 2], [-69, -260], [134, -227], [-41, -278], [-65, -27]], [[55338, 76294], [-52, -53], [-90, -138], [-41, -325]], [[55155, 75778], [-246, 224], [-105, 247], [-106, 130], [-127, 221], [-61, 183], [-136, 277], [59, 245], [99, -136], [60, 123], [130, 13], [239, -98], [192, 8], [126, -131]], [[56523, 82432], [268, -4], [302, 223], [64, 333], [228, 190], [-26, 264]], [[57359, 83438], [169, 100], [298, 228]], [[57826, 83766], [293, -149], [39, -146], [146, 70], [272, -141], [27, -277], [-60, -159], [174, -387], [113, -108], [-16, -107], [187, -104], [80, -157], [-108, -129], [-224, 20], [-54, -55], [66, -196], [68, -379]], [[58829, 81362], [-239, -35], [-85, -129], [-18, -298], [-111, 57], [-250, -28], [-73, 138], [-104, -103], [-105, 86], [-218, 12], [-310, 141], [-281, 47], [-215, -14], [-152, -160], [-133, -23]], [[56535, 81053], [-6, 263], [-85, 274], [166, 121], [2, 235], [-77, 225], [-12, 261]], [[25238, 61101], [-2, 87], [33, 27], [51, -70], [99, 357], [53, 8]], [[25472, 61510], [1, -87], [53, -3], [-5, -160], [-45, -256], [24, -91], [-29, -212], [18, -56], [-32, -299], [-55, -156], [-50, -19], [-55, -205]], [[25297, 59966], [-83, 0], [22, 667], [2, 468]], [[31359, 37147], [-200, -81], [-109, 814], [-150, 663], [88, 572], [-146, 250], [-37, 426], [-136, 402]], [[30669, 40193], [175, 638], [-119, 496], [63, 199], [-49, 219], [108, 295], [6, 503], [13, 415], [60, 200], [-240, 951]], [[30686, 44109], [206, -50], [143, 13], [62, 179], [243, 239], [147, 222], [363, 100], [-29, -443], [34, -227], [-23, -396], [302, -529], [311, -98], [109, -220], [188, -117], [115, -172], [175, 6], [161, -175], [12, -342], [55, -172], [3, -255], [-81, -10], [107, -688], [533, -24], [-41, -342], [30, -233], [151, -166], [66, -367], [-49, -465], [-77, -259], [27, -337], [-87, -122]], [[33842, 38659], [-4, 182], [-259, 302], [-258, 9], [-484, -172], [-133, -520], [-7, -318], [-110, -708]], [[34826, 35372], [54, 341], [38, 350], [0, 325], [-100, 107], [-104, -96], [-103, 26], [-33, 228], [-26, 541], [-52, 177], [-187, 160], [-114, -116], [-293, 113], [18, 802], [-82, 329]], [[30686, 44109], [-157, -102], [-126, 68], [18, 898], [-228, -348], [-245, 15], [-105, 315], [-184, 34], [59, 254], [-155, 359], [-115, 532], [73, 108], [0, 250], [168, 171], [-28, 319], [71, 206], [20, 275], [318, 402], [227, 114], [37, 89], [251, -28]], [[30585, 48040], [125, 1620], [6, 256], [-43, 339], [-123, 215], [1, 430], [156, 97], [56, -61], [9, 226], [-162, 61], [-4, 370], [541, -13], [92, 203], [77, -187], [55, -349], [52, 73]], [[31423, 51320], [153, -312], [216, 38], [54, 181], [206, 138], [115, 97], [32, 250], [198, 168], [-15, 124], [-235, 51], [-39, 372], [12, 396], [-125, 153], [52, 55], [206, -76], [221, -148], [80, 140], [200, 92], [310, 221], [102, 225], [-37, 167]], [[33129, 53652], [145, 26], [64, -136], [-36, -259], [96, -90], [63, -274], [-77, -209], [-44, -502], [71, -299], [20, -274], [171, -277], [137, -29], [30, 116], [88, 25], [126, 104], [90, 157], [154, -50], [67, 21]], [[34294, 51702], [151, -48], [25, 120], [-46, 118], [28, 171], [112, -53], [131, 61], [159, -125]], [[34854, 51946], [121, -122], [86, 160], [62, -25], [38, -166], [133, 42], [107, 224], [85, 436], [164, 540]], [[35650, 53035], [95, 28], [69, -327], [155, -1033], [149, -97], [7, -408], [-208, -487], [86, -178], [491, -92], [10, -593], [211, 388], [349, -212], [462, -361], [135, -346], [-45, -327], [323, 182], [540, -313], [415, 23], [411, -489], [355, -662], [214, -170], [237, -24], [101, -186], [94, -752], [46, -358], [-110, -977], [-142, -385], [-391, -822], [-177, -668], [-206, -513], [-69, -11], [-78, -435], [20, -1107], [-77, -910], [-30, -390], [-88, -233], [-49, -790], [-282, -771], [-47, -610], [-225, -256], [-65, -355], [-302, 2], [-437, -227], [-195, -263], [-311, -173], [-327, -470], [-235, -586], [-41, -441], [46, -326], [-51, -597], [-63, -289], [-195, -325], [-308, -1040], [-244, -468], [-189, -277], [-127, -562], [-183, -337]], [[35174, 30629], [-77, 334], [122, 280], [-160, 402], [-218, 327], [-286, 379], [-103, -18], [-279, 457], [-180, -63]], [[81723, 53254], [110, 221], [236, 323]], [[82069, 53798], [-13, -291], [-16, -377], [-133, 19], [-58, -202], [-126, 307]], [[75471, 66988], [113, -189], [-20, -363], [-227, -17], [-234, 39], [-175, -92], [-252, 224], [-6, 119]], [[74670, 66709], [184, 439], [150, 150], [198, -137], [147, -14], [122, -159]], [[58175, 37528], [-393, -435], [-249, -442], [-93, -393], [-83, -222], [-152, -47], [-48, -283], [-28, -184], [-178, -138], [-226, 29], [-133, 166], [-117, 71], [-135, -137], [-68, -283], [-132, -177], [-139, -264], [-199, -60], [-62, 207], [26, 360], [-165, 562], [-75, 88]], [[55526, 35946], [0, 1725], [274, 20], [8, 2105], [207, 19], [428, 207], [106, -243], [177, 231], [85, 2], [156, 133]], [[56967, 40145], [50, -44]], [[57017, 40101], [107, -473], [56, -105], [87, -342], [315, -649], [119, -64], [0, -208], [82, -375], [215, -90], [177, -267]], [[54244, 54965], [229, 44], [52, 152], [46, -11], [69, -134], [350, 226], [118, 230], [145, 207], [-28, 208], [78, 54], [269, -36], [261, 273], [201, 645], [141, 239], [176, 101]], [[56351, 57163], [31, -253], [160, -369], [1, -241], [-45, -246], [18, -184], [96, -170]], [[56612, 55700], [212, -258]], [[56824, 55442], [152, -239], [2, -192], [187, -308], [116, -255], [70, -355], [208, -234], [44, -187]], [[57603, 53672], [-91, -63], [-178, 14], [-209, 62], [-104, -51], [-41, -143], [-90, -18], [-110, 125], [-309, -295], [-127, 60], [-38, -46], [-83, -357], [-207, 115], [-203, 59], [-177, 218], [-229, 200], [-149, -190], [-108, -300], [-25, -412]], [[55125, 52650], [-178, 33], [-188, 99], [-166, -313], [-146, -550]], [[54447, 51919], [-29, 172], [-12, 269], [-127, 190], [-103, 305], [-23, 212], [-132, 309], [23, 176], [-28, 249], [21, 458], [67, 107], [140, 599]], [[32315, 78082], [202, -79], [257, 16], [-137, -242], [-102, -38], [-353, 250], [-69, 198], [105, 183], [97, -288]], [[32831, 79592], [-135, -11], [-360, 186], [-258, 279], [96, 49], [365, -148], [284, -247], [8, -108]], [[15692, 79240], [-140, -82], [-456, 269], [-84, 209], [-248, 207], [-50, 168], [-286, 107], [-107, 321], [24, 137], [291, -129], [171, -89], [261, -63], [94, -204], [138, -280], [277, -244], [115, -327]], [[34407, 80527], [-184, -517], [181, 199], [187, -126], [-98, -206], [247, -162], [128, 144], [277, -182], [-86, -433], [194, 101], [36, -313], [86, -367], [-117, -520], [-125, -22], [-183, 111], [60, 484], [-77, 75], [-322, -513], [-166, 21], [196, 277], [-267, 144], [-298, -35], [-539, 18], [-43, 175], [173, 208], [-121, 160], [234, 356], [287, 941], [172, 336], [241, 204], [129, -26], [-54, -160], [-148, -372]], [[13005, 82584], [131, -76], [267, 47], [-84, -671], [242, -475], [-111, 1], [-167, 270], [-103, 272], [-140, 184], [-51, 260], [16, 188]], [[27981, 87304], [-108, -310], [-123, 50], [-73, 176], [13, 41], [107, 177], [114, -13], [70, -121]], [[27250, 87631], [-325, -326], [-196, 13], [-61, 160], [207, 273], [381, -6], [-6, -114]], [[26344, 89371], [51, -259], [143, 91], [161, -155], [304, -203], [318, -184], [25, -281], [204, 46], [199, -196], [-247, -186], [-432, 142], [-156, 266], [-275, -314], [-396, -306], [-95, 346], [-377, -57], [242, 292], [35, 465], [95, 542], [201, -49]], [[28926, 90253], [-312, -30], [-69, 289], [118, 331], [255, 82], [217, -163], [3, -253], [-32, -82], [-180, -174]], [[23431, 91410], [-173, -207], [-374, 179], [-226, -65], [-380, 266], [245, 183], [194, 256], [295, -168], [166, -106], [84, -112], [169, -226]], [[31350, 77248], [-181, 334], [0, 805], [-123, 171], [-187, -100], [-92, 155], [-212, -446], [-84, -460], [-99, -269], [-118, -91], [-89, -30], [-28, -146], [-512, 0], [-422, -4], [-125, -109], [-294, -425], [-34, -46], [-89, -231], [-255, 1], [-273, -3], [-125, -93], [44, -116], [25, -181], [-5, -60], [-363, -293], [-286, -93], [-323, -316], [-70, 0], [-94, 93], [-31, 85], [6, 61], [61, 207], [131, 325], [81, 349], [-56, 514], [-59, 536], [-290, 277], [35, 105], [-41, 73], [-76, 0], [-56, 93], [-14, 140], [-54, -61], [-75, 18], [17, 59], [-65, 58], [-27, 155], [-216, 189], [-224, 197], [-272, 229], [-261, 214], [-248, -167], [-91, -6], [-342, 154], [-225, -77], [-269, 183], [-284, 94], [-194, 36], [-86, 100], [-49, 325], [-94, -3], [-1, -227], [-575, 0], [-951, 0], [-944, 0], [-833, 0], [-834, 0], [-819, 0], [-847, 0], [-273, 0], [-825, 0], [-788, 0]], [[15878, 79530], [-38, 1], [-537, 581], [-199, 255], [-503, 244], [-155, 523], [40, 363], [-356, 252], [-48, 476], [-336, 429], [-6, 304]], [[13740, 82958], [154, 285], [-7, 373], [-473, 376], [-284, 674], [-173, 424], [-255, 266], [-187, 242], [-147, 306], [-279, -192], [-270, -330], [-247, 388], [-194, 259], [-271, 164], [-273, 17], [1, 3364], [2, 2193]], [[10837, 91767], [518, -142], [438, -285], [289, -54], [244, 247], [336, 184], [413, -72], [416, 259], [455, 148], [191, -245], [207, 138], [62, 278], [192, -63], [470, -530], [369, 401], [38, -449], [341, 97], [105, 173], [337, -34], [424, -248], [650, -217], [383, -100], [272, 38], [374, -300], [-390, -293], [502, -127], [750, 70], [236, 103], [296, -354], [302, 299], [-283, 251], [179, 202], [338, 27], [223, 59], [224, -141], [279, -321], [310, 47], [491, -266], [431, 94], [405, -14], [-32, 367], [247, 103], [431, -200], [-2, -559], [177, 471], [223, -16], [126, 594], [-298, 364], [-324, 239], [22, 653], [329, 429], [366, -95], [281, -261], [378, -666], [-247, -290], [517, -120], [-1, -604], [371, 463], [332, -380], [-83, -438], [269, -399], [290, 427], [202, 510], [16, 649], [394, -46], [411, -87], [373, -293], [17, -293], [-207, -315], [196, -316], [-36, -288], [-544, -413], [-386, -91], [-287, 178], [-83, -297], [-268, -498], [-81, -259], [-322, -399], [-397, -39], [-220, -250], [-18, -384], [-323, -74], [-340, -479], [-301, -665], [-108, -466], [-16, -686], [409, -99], [125, -553], [130, -448], [388, 117], [517, -256], [277, -225], [199, -279], [348, -163], [294, -248], [459, -34], [302, -58], [-45, -511], [86, -594], [201, -661], [414, -561], [214, 192], [150, 607], [-145, 934], [-196, 311], [445, 276], [314, 415], [154, 411], [-23, 395], [-188, 502], [-338, 445], [328, 619], [-121, 535], [-93, 922], [194, 137], [476, -161], [286, -57], [230, 155], [258, -200], [342, -343], [85, -229], [495, -45], [-8, -496], [92, -747], [254, -92], [201, -348], [402, 328], [266, 652], [184, 274], [216, -527], [362, -754], [307, -709], [-112, -371], [370, -333], [250, -338], [442, -152], [179, -189], [110, -500], [216, -78], [112, -223], [20, -664], [-202, -222], [-199, -207], [-458, -210], [-349, -486], [-470, -96], [-594, 125], [-417, 4], [-287, -41], [-233, -424], [-354, -262], [-401, -782], [-320, -545], [236, 97], [446, 776], [583, 493], [415, 58], [246, -289], [-262, -397], [88, -637], [91, -446], [361, -295], [459, 86], [278, 664], [19, -429], [180, -214], [-344, -387], [-615, -351], [-276, -239], [-310, -426], [-211, 44], [-11, 500], [483, 488], [-445, -19], [-309, -72]], [[18287, 93781], [-139, -277], [618, 179], [386, -298], [314, 302], [254, -194], [227, -580], [140, 244], [-197, 606], [244, 86], [276, -94], [311, -239], [175, -575], [86, -417], [466, -293], [502, -279], [-31, -260], [-456, -48], [178, -227], [-94, -217], [-503, 93], [-478, 160], [-322, -36], [-522, -201], [-704, -88], [-494, -56], [-151, 279], [-379, 161], [-246, -66], [-343, 468], [185, 62], [429, 101], [392, -26], [362, 103], [-537, 138], [-594, -47], [-394, 12], [-146, 217], [644, 237], [-428, -9], [-485, 156], [233, 443], [193, 235], [744, 359], [284, -114]], [[20972, 93958], [-244, -390], [-434, 413], [95, 83], [372, 24], [211, -130]], [[28794, 93770], [25, -163], [-296, 17], [-299, 13], [-304, -80], [-80, 36], [-306, 313], [12, 213], [133, 39], [636, -63], [479, -325]], [[25955, 93803], [219, -369], [256, 477], [704, 242], [477, -611], [-42, -387], [550, 172], [263, 235], [616, -299], [383, -282], [36, -258], [515, 134], [290, -376], [670, -234], [242, -238], [263, -553], [-510, -275], [654, -386], [441, -130], [400, -543], [437, -39], [-87, -414], [-487, -687], [-342, 253], [-437, 568], [-359, -74], [-35, -338], [292, -344], [377, -272], [114, -157], [181, -584], [-96, -425], [-350, 160], [-697, 473], [393, -509], [289, -357], [45, -206], [-753, 236], [-596, 343], [-337, 287], [97, 167], [-414, 304], [-405, 286], [5, -171], [-803, -94], [-235, 203], [183, 435], [522, 10], [571, 76], [-92, 211], [96, 294], [360, 576], [-77, 261], [-107, 203], [-425, 286], [-563, 201], [178, 150], [-294, 367], [-245, 34], [-219, 201], [-149, -175], [-503, -76], [-1011, 132], [-588, 174], [-450, 89], [-231, 207], [290, 270], [-394, 2], [-88, 599], [213, 528], [286, 241], [717, 158], [-204, -382]], [[22123, 94208], [331, -124], [496, 75], [72, -172], [-259, -283], [420, -254], [-50, -532], [-455, -229], [-268, 50], [-192, 225], [-690, 456], [5, 189], [567, -73], [-306, 386], [329, 286]], [[24112, 93575], [-298, -442], [-317, 22], [-173, 519], [4, 294], [145, 251], [276, 161], [579, -20], [530, -144], [-415, -526], [-331, -115]], [[16539, 92755], [-731, -285], [-147, 259], [-641, 312], [119, 250], [192, 432], [241, 388], [-272, 362], [939, 93], [397, -123], [709, -33], [270, -171], [298, -249], [-349, -149], [-681, -415], [-344, -414], [0, -257]], [[23996, 94879], [-151, -229], [-403, 44], [-337, 155], [148, 266], [399, 159], [243, -208], [101, -187]], [[22639, 95907], [212, -273], [9, -303], [-127, -440], [-458, -60], [-298, 94], [5, 345], [-455, -46], [-18, 457], [299, -18], [419, 201], [390, -34], [22, 77]], [[19941, 95601], [109, -210], [247, 99], [291, -26], [49, -289], [-169, -281], [-940, -91], [-701, -256], [-423, -14], [-35, 193], [577, 261], [-1255, -70], [-389, 106], [379, 577], [262, 165], [782, -199], [493, -350], [485, -45], [-397, 565], [255, 215], [286, -68], [94, -282]], [[23699, 96131], [308, -190], [547, 1], [240, -194], [-64, -222], [319, -134], [177, -140], [374, -26], [406, -50], [441, 128], [566, 51], [451, -42], [298, -223], [62, -244], [-174, -157], [-414, -127], [-355, 72], [-797, -91], [-570, -11], [-449, 73], [-738, 190], [-96, 325], [-34, 293], [-279, 258], [-574, 72], [-322, 183], [104, 242], [573, -37]], [[17722, 96454], [-38, -454], [-214, -205], [-259, -29], [-517, -252], [-444, -91], [-377, 128], [472, 442], [570, 383], [426, -9], [381, 87]], [[23933, 96380], [-126, -17], [-521, 38], [-74, 165], [559, -9], [195, -109], [-33, -68]], [[19392, 96485], [-518, -170], [-411, 191], [224, 188], [406, 60], [392, -92], [-93, -177]], [[19538, 97019], [-339, -115], [-461, 1], [5, 84], [285, 177], [149, -27], [361, -120]], [[23380, 96697], [-411, -122], [-226, 138], [-119, 221], [-22, 245], [360, -24], [162, -39], [332, -205], [-76, -214]], [[22205, 96856], [108, -247], [-453, 66], [-457, 192], [-619, 21], [268, 176], [-335, 142], [-21, 227], [546, -81], [751, -215], [212, -281]], [[25828, 97644], [334, -190], [-381, -176], [-513, -445], [-492, -42], [-575, 76], [-299, 240], [4, 215], [220, 157], [-508, -4], [-306, 196], [-176, 268], [193, 262], [192, 180], [285, 42], [-122, 135], [646, 30], [355, -315], [468, -127], [455, -112], [220, -390]], [[30972, 99681], [742, -47], [597, -75], [508, -161], [-12, -157], [-678, -257], [-672, -119], [-251, -133], [605, 3], [-656, -358], [-452, -167], [-476, -483], [-573, -98], [-177, -120], [-841, -64], [383, -74], [-192, -105], [230, -292], [-264, -202], [-429, -167], [-132, -232], [-388, -176], [39, -134], [475, 23], [6, -144], [-742, -355], [-726, 163], [-816, -91], [-414, 71], [-525, 31], [-35, 284], [514, 133], [-137, 427], [170, 41], [742, -255], [-379, 379], [-450, 113], [225, 229], [492, 141], [79, 206], [-392, 231], [-118, 304], [759, -26], [220, -64], [433, 216], [-625, 68], [-972, -38], [-491, 201], [-232, 239], [-324, 173], [-61, 202], [413, 112], [324, 19], [545, 96], [409, 220], [344, -30], [300, -166], [211, 319], [367, 95], [498, 65], [849, 24], [148, -63], [802, 100], [601, -38], [602, -37]], [[52900, 78285], [-22, -242], [-122, -100], [-206, 75], [-60, -239], [-132, -19], [-48, 94], [-156, -200], [-134, -28], [-120, 126]], [[51900, 77752], [-95, 259], [-133, -92], [5, 267], [203, 332], [-9, 150], [126, -54], [77, 101]], [[52074, 78715], [236, -4], [57, 128], [298, -181]], [[31400, 18145], [-92, -239], [-238, -183], [-137, 19], [-164, 48], [-202, 177], [-291, 86], [-350, 330], [-283, 317], [-383, 662], [229, -124], [390, -395], [369, -212], [143, 271], [90, 405], [256, 244], [198, -70]], [[30952, 19680], [-247, 4], [-134, -145], [-250, -213], [-45, -552], [-118, -14], [-313, 192], [-318, 412], [-346, 338], [-87, 374], [79, 346], [-140, 393], [-36, 1007], [119, 568], [293, 457], [-422, 172], [265, 522], [94, 982], [309, -208], [145, 1224], [-186, 157], [-87, -738], [-175, 83], [87, 845], [95, 1095], [127, 404], [-80, 576], [-22, 666], [117, 19], [170, 954], [192, 945], [118, 881], [-64, 885], [83, 487], [-34, 730], [163, 721], [50, 1143], [89, 1227], [87, 1321], [-20, 967], [-58, 832]], [[30452, 39739], [143, 151], [74, 303]], [[80649, 61615], [-240, -284], [-228, 183], [-8, 509], [137, 267], [304, 166], [159, -14], [62, -226], [-122, -260], [-64, -341]], [[86288, 75628], [-179, 348], [-111, -331], [-429, -254], [44, -312], [-241, 22], [-131, 185], [-191, -419], [-306, -318], [-227, -379]], [[84517, 74170], [-388, -171], [-204, -277], [-300, -161], [148, 274], [-58, 230], [220, 397], [-147, 310], [-242, -209], [-314, -411], [-171, -381], [-272, -29], [-142, -275], [147, -400], [227, -97], [9, -265], [220, -173], [311, 422], [247, -230], [179, -15], [45, -310], [-393, -165], [-130, -319], [-270, -296], [-142, -414], [299, -325], [109, -581], [169, -541], [189, -454], [-5, -439], [-174, -161], [66, -315], [164, -184], [-43, -481], [-71, -468], [-155, -53], [-203, -640], [-225, -775], [-258, -705], [-382, -545], [-386, -498], [-313, -68], [-170, -262], [-96, 192], [-157, -294], [-388, -296], [-294, -90], [-95, -624], [-154, -35], [-73, 429], [66, 228], [-373, 189], [-131, -96]], [[80013, 63313], [-280, 154], [-132, 240], [44, 340], [-254, 108], [-134, 222], [-236, -315], [-271, -68], [-221, 3], [-149, -145]], [[78380, 63852], [-144, -86], [42, -676], [-148, 16], [-25, 139]], [[78105, 63245], [-9, 244], [-203, -172], [-121, 109], [-206, 222], [81, 490], [-176, 115], [-66, 544], [-293, -98], [33, 701], [263, 493], [11, 487], [-8, 452], [-121, 141], [-93, 348], [-162, -44]], [[77035, 67277], [-300, 89], [94, 248], [-130, 367], [-198, -249], [-233, 145], [-321, -376], [-252, -439], [-224, -74]], [[74670, 66709], [-23, 465], [-170, -124]], [[74477, 67050], [-324, 57], [-314, 136], [-225, 259], [-216, 117], [-93, 284], [-157, 84], [-280, 385], [-223, 182], [-115, -141]], [[72530, 68413], [-386, 413], [-273, 374], [-78, 651], [200, -79], [9, 301], [-111, 303], [28, 482], [-298, 692]], [[71621, 71550], [-457, 239], [-82, 454], [-205, 276]], [[70827, 72688], [-42, 337], [10, 230], [-169, 134], [-91, -59], [-70, 546]], [[70465, 73876], [79, 136], [-39, 138], [266, 279], [192, 116], [294, -80], [105, 378], [356, 70], [99, 234], [438, 320], [39, 134]], [[72294, 75601], [-22, 337], [190, 154], [-250, 1026], [550, 236], [143, 131], [200, 1058], [551, -194], [155, 267], [13, 592], [230, 56], [212, 393]], [[74266, 79657], [109, 49]], [[74375, 79706], [73, -413], [233, -313], [396, -222], [192, -476], [-107, -690], [100, -256], [330, -101], [374, -83], [336, -368], [171, -66], [127, -544], [163, -351], [306, 14], [574, -133], [369, 82], [274, -88], [411, -359], [336, 1], [123, -184], [324, 318], [448, 205], [417, 22], [324, 208], [200, 316], [194, 199], [-45, 195], [-89, 227], [146, 381], [156, -53], [286, -120], [277, 313], [423, 229], [204, 391], [195, 168], [404, 78], [219, -66], [30, 210], [-251, 413], [-223, 189], [-214, -219], [-274, 92], [-157, -74], [-72, 241], [197, 590], [135, 446]], [[82410, 80055], [333, -223], [392, 373], [-3, 260], [251, 627], [155, 189], [-4, 326], [-152, 141], [229, 294], [345, 106], [369, 16], [415, -176], [244, -217], [172, -596], [104, -254], [97, -363], [103, -579], [483, -189], [329, -420], [112, -555], [423, -1], [240, 233], [459, 175], [-146, -532], [-107, -216], [-96, -647], [-186, -575], [-338, 104], [-238, -208], [73, -506], [-40, -698], [-142, -16], [2, -300]], [[49206, 53531], [-126, -7], [-194, 116], [-178, -7], [-329, -103], [-193, -170], [-275, -217], [-54, 15]], [[47857, 53158], [22, 487], [26, 74], [-8, 233], [-118, 247], [-88, 40], [-81, 162], [60, 262], [-28, 286], [13, 172]], [[47655, 55121], [44, 0], [17, 258], [-22, 114], [27, 82], [103, 71], [-69, 473], [-64, 245], [23, 200], [55, 46]], [[47769, 56610], [36, 54], [77, -89], [215, -5], [51, 172], [48, -11], [80, 67], [43, -253], [65, 74], [114, 88]], [[49214, 56277], [74, -841], [-117, -496], [-73, -667], [121, -509], [-13, -233]], [[53632, 51919], [-35, 32], [-164, -76], [-169, 79], [-132, -38]], [[53132, 51916], [-452, 13]], [[52680, 51929], [40, 466], [-108, 391], [-127, 100], [-56, 265], [-72, 85], [4, 163]], [[52361, 53399], [71, 418], [132, 570], [81, 6], [165, 345], [105, 10], [156, -243], [191, 199], [26, 246], [63, 238], [43, 299], [148, 243], [56, 414], [59, 132], [39, 307], [74, 377], [234, 457], [14, 196], [31, 107], [-110, 235]], [[53939, 57955], [9, 188], [78, 34]], [[54026, 58177], [111, -378], [18, -392], [-10, -393], [151, -537], [-155, 6], [-78, -42], [-127, 60], [-60, -279], [164, -345], [121, -100], [39, -245], [87, -407], [-43, -160]], [[54447, 51919], [-20, -319], [-220, 140], [-225, 156], [-350, 23]], [[58564, 52653], [-16, -691], [111, -80], [-89, -210], [-107, -157], [-106, -308], [-59, -274], [-15, -475], [-65, -225], [-2, -446]], [[58216, 49787], [-80, -165], [-10, -351], [-38, -46], [-26, -323]], [[58149, 47921], [50, -544], [-27, -307]], [[58172, 47070], [55, -343], [161, -330]], [[58388, 46397], [150, -745]], [[58538, 45652], [-109, 60], [-373, -99], [-75, -71], [-79, -377], [62, -261], [-49, -699], [-34, -593], [75, -105], [194, -230], [76, 107], [23, -637], [-212, 5], [-114, 325], [-103, 252], [-213, 82], [-62, 310], [-170, -187], [-222, 83], [-93, 268], [-176, 55], [-131, -15], [-15, 184], [-96, 15]], [[53422, 46976], [-39, 183]], [[53609, 47755], [73, -60], [95, 226], [152, -6], [17, -167], [104, -105], [164, 370], [161, 289], [71, 189], [-10, 486], [121, 574], [127, 304], [183, 285], [32, 189], [7, 216], [45, 205], [-14, 335], [34, 524], [55, 368], [83, 316], [16, 357]], [[57603, 53672], [169, -488], [124, -71], [75, 99], [128, -39], [155, 125], [66, -252], [244, -393]], [[53309, 47603], [-228, 626]], [[53081, 48229], [212, 326], [-105, 391], [95, 148], [187, 73], [23, 261], [148, -283], [245, -25], [85, 279], [36, 393], [-31, 461], [-131, 350], [120, 684], [-69, 117], [-207, -48], [-78, 305], [21, 258]], [[29063, 50490], [-119, 140], [-137, 195], [-79, -94], [-235, 82], [-68, 255], [-52, -10], [-278, 338]], [[28095, 51396], [-37, 183], [103, 44], [-12, 296], [65, 214], [138, 40], [117, 371], [106, 310], [-102, 141], [52, 343], [-62, 540], [59, 155], [-44, 500], [-112, 315]], [[28366, 54848], [36, 287], [89, -43], [52, 176], [-64, 348], [34, 86]], [[28513, 55702], [143, -18], [209, 412], [114, 63], [3, 195], [51, 500], [159, 274], [175, 11], [22, 123], [218, -49], [218, 298], [109, 132], [134, 285], [98, -36], [73, -156], [-54, -199]], [[30185, 57537], [-178, -99], [-71, -295], [-107, -169], [-81, -220], [-34, -422], [-77, -345], [144, -40], [35, -271], [62, -130], [21, -238], [-33, -219], [10, -123], [69, -49], [66, -207], [357, 57], [161, -75], [196, -508], [112, 63], [200, -32], [158, 68], [99, -102], [-50, -318], [-62, -199], [-22, -423], [56, -393], [79, -175], [9, -133], [-140, -294], [100, -130], [74, -207], [85, -589]], [[30585, 48040], [-139, 314], [-83, 14], [179, 602], [-213, 276], [-166, -51], [-101, 103], [-153, -157], [-207, 74], [-163, 620], [-129, 152], [-89, 279], [-184, 280], [-74, -56]], [[26954, 55439], [-151, 131], [-56, 124], [32, 103], [-11, 130], [-77, 142], [-109, 116], [-95, 76], [-19, 173], [-73, 105], [18, -172], [-55, -141], [-64, 164], [-89, 58], [-38, 120], [2, 179], [36, 187], [-78, 83], [64, 114]], [[26191, 57131], [42, 76], [183, -156], [63, 77], [89, -50], [46, -121], [82, -40], [66, 126]], [[26762, 57043], [70, -321], [108, -238], [130, -252]], [[27070, 56232], [-107, -53], [1, -238], [58, -88], [-41, -70], [10, -107], [-23, -120], [-14, -117]], [[27147, 64280], [240, -42], [219, -7], [261, -201], [110, -216], [260, 66], [98, -138], [235, -366], [173, -267], [92, 8], [165, -120], [-20, -167], [205, -24], [210, -242], [-33, -138], [-185, -75], [-187, -29], [-191, 46], [-398, -57], [186, 329], [-113, 154], [-179, 39], [-96, 171], [-66, 336], [-157, -23], [-259, 159], [-83, 124], [-362, 91], [-97, 115], [104, 148], [-273, 30], [-199, -307], [-115, -8], [-40, -144], [-138, -65], [-118, 56], [146, 183], [60, 213], [126, 131], [142, 116], [210, 56], [67, 65]], [[59092, 71341], [19, 3], [40, 143], [200, -8], [253, 176], [-188, -251], [21, -111]], [[59437, 71293], [-30, 21], [-53, -45], [-42, 12], [-14, -22], [-5, 59], [-20, 37], [-54, 6], [-75, -51], [-52, 31]], [[59437, 71293], [8, -48], [-285, -240], [-136, 77], [-64, 237], [132, 22]], [[53776, 79457], [-157, 254], [-141, 142], [-30, 249], [-49, 176], [202, 129], [103, 147], [200, 114], [70, 113], [73, -68], [124, 62]], [[54171, 80775], [132, -191], [207, -51], [-17, -163], [151, -122], [41, 153], [191, -66], [26, -185], [207, -36], [127, -291]], [[55236, 79823], [-82, -1], [-43, -106], [-64, -26], [-18, -134], [-54, -28], [-7, -55], [-95, -61], [-123, 10], [-39, -130]], [[52756, 83065], [4, -228], [281, -138], [-3, -210], [283, 111], [156, 162], [313, -233], [132, -189]], [[53922, 82340], [64, -300], [-77, -158], [101, -210], [69, -316], [-22, -204], [114, -377]], [[52074, 78715], [35, 421], [140, 404], [-400, 109], [-131, 155]], [[51718, 79804], [16, 259], [-56, 133]], [[51710, 80596], [-47, 619], [167, 0], [70, 222], [69, 541], [-51, 200]], [[51918, 82178], [54, 125], [232, 32], [52, -130], [188, 291], [-63, 222], [-13, 335]], [[52368, 83053], [210, -78], [178, 90]], [[61966, 58083], [66, -183], [-9, -245], [-158, -142], [119, -161]], [[61984, 57352], [-102, -317]], [[61882, 57035], [-62, 106], [-67, -42], [-155, 10], [-4, 180], [-22, 163], [94, 277], [98, 261]], [[61764, 57990], [119, -51], [83, 144]], [[53524, 83435], [-166, -478], [-291, 333], [-39, 246], [408, 195], [88, -296]], [[52368, 83053], [-113, 328], [-8, 604], [46, 159], [80, 177], [244, 37], [98, 163], [223, 167], [-9, -304], [-82, -192], [33, -166], [151, -89], [-68, -223], [-83, 64], [-200, -425], [76, -288]], [[30080, 62227], [34, 101], [217, -3], [165, -152], [73, 15], [50, -209], [152, 11], [-9, -176], [124, -21], [136, -217], [-103, -240], [-132, 128], [-127, -25], [-92, 28], [-50, -107], [-106, -37], [-43, 144], [-92, -85], [-111, -405], [-71, 94], [-14, 170]], [[30081, 61241], [5, 161], [-71, 177], [68, 99], [21, 228], [-24, 321]], [[53333, 64447], [-952, -1126], [-804, -1161], [-392, -263]], [[51185, 61897], [-308, -58], [-3, 376], [-129, 96], [-173, 169], [-66, 277], [-937, 1289], [-937, 1289]], [[48632, 65335], [-1045, 1431]], [[47587, 66766], [6, 114], [-1, 40]], [[47592, 66920], [-2, 700], [449, 436], [277, 90], [227, 159], [107, 295], [324, 234], [12, 438], [161, 51], [126, 219], [363, 99], [51, 230], [-73, 125], [-96, 624], [-17, 359], [-104, 379]], [[49397, 71358], [267, 323], [300, 102], [175, 244], [268, 180], [471, 105], [459, 48], [140, -87], [262, 232], [297, 5], [113, -137], [190, 35]], [[52339, 72408], [-57, -303], [44, -563], [-65, -487], [-171, -330], [24, -445], [227, -352], [3, -143], [171, -238], [118, -1061]], [[52633, 68486], [90, -522], [15, -274], [-49, -482], [21, -270], [-36, -323], [24, -371], [-110, -247], [164, -431], [11, -253], [99, -330], [130, 109], [219, -275], [122, -370]], [[27693, 48568], [148, 442], [-60, 258], [-106, -275], [-166, 259], [56, 167], [-47, 536], [97, 89], [52, 368], [105, 381], [-20, 241], [153, 126], [190, 236]], [[29063, 50490], [38, -449], [-86, -384], [-303, -619], [-334, -233], [-170, -514], [-53, -398], [-157, -243], [-116, 298], [-113, 64], [-114, -47], [-8, 216], [79, 141], [-33, 246]], [[59700, 68010], [-78, -238], [-60, -446], [-75, -308], [-65, -103], [-93, 191], [-125, 263], [-198, 847], [-29, -53], [115, -624], [171, -594], [210, -920], [102, -321], [90, -334], [249, -654], [-55, -103], [9, -384], [323, -530], [49, -121]], [[60240, 63578], [-1102, 0], [-1077, 0], [-1117, 0]], [[56944, 63578], [0, 2175], [0, 2101], [-83, 476], [71, 365], [-43, 253], [101, 283]], [[56990, 69231], [369, 10], [268, -156], [275, -175], [129, -92], [214, 188], [114, 169], [245, 49], [198, -75], [75, -293], [65, 193], [222, -140], [217, -33], [137, 149]], [[59518, 69025], [182, -1015]], [[61764, 57990], [-95, 191], [-114, 346], [-124, 190], [-71, 204], [-242, 237], [-191, 7], [-67, 124], [-163, -139], [-168, 268], [-87, -441], [-323, 124]], [[60119, 59101], [-30, 236], [120, 868], [27, 393], [88, 181], [204, 97], [141, 337]], [[60669, 61213], [161, -684], [77, -542], [152, -288], [379, -558], [154, -336], [151, -341], [87, -203], [136, -178]], [[47490, 75324], [14, 420], [-114, 257], [393, 426], [340, -106], [373, 3], [296, -101], [230, 31], [449, -19]], [[49471, 76235], [111, -230], [511, -268], [101, 127], [313, -267], [322, 77]], [[50829, 75674], [15, -344], [-263, -393], [-356, -125], [-25, -199], [-171, -327], [-107, -481], [108, -338], [-160, -263], [-60, -384], [-210, -118], [-197, -454], [-352, -9], [-265, 11], [-174, -209], [-106, -223], [-136, 49], [-103, 199], [-79, 340], [-259, 92]], [[47929, 72498], [-23, 195], [103, 222], [38, 161], [-96, 175], [77, 388], [-111, 355], [120, 48], [11, 280], [45, 86], [3, 461], [129, 160], [-78, 296], [-162, 21], [-47, -75], [-164, 0], [-70, 289], [-113, -86], [-101, -150]], [[56753, 84725], [32, 349], [-102, -75], [-176, 210], [-24, 340], [351, 164], [350, 86], [301, -97], [287, 17]], [[57772, 85719], [42, -103], [-198, -341], [83, -551], [-120, -187]], [[57579, 84537], [-229, 1], [-239, 219], [-121, 73], [-237, -105]], [[61882, 57035], [-61, -209], [103, -325], [102, -285], [106, -210], [909, -702], [233, 4]], [[63274, 55308], [-785, -1773], [-362, -26], [-247, -417], [-178, -11], [-76, -186]], [[61626, 52895], [-190, 0], [-112, 200], [-254, -247], [-82, -247], [-185, 47], [-62, 68], [-65, -16], [-87, 6], [-352, 502], [-193, 0], [-95, 194], [0, 332], [-145, 99]], [[59804, 53833], [-164, 643], [-127, 137], [-48, 236], [-141, 288], [-171, 42], [95, 337], [147, 14], [42, 181]], [[59437, 55711], [-4, 531]], [[59433, 56242], [82, 618], [132, 166], [28, 241], [119, 451], [168, 293], [112, 582], [45, 508]], [[57942, 91385], [-41, -414], [425, -394], [-256, -445], [323, -673], [-187, -506], [250, -440], [-113, -385], [411, -405], [-105, -301], [-258, -341], [-594, -755]], [[57797, 86326], [-504, -47], [-489, -216], [-452, -125], [-161, 323], [-269, 193], [62, 582], [-135, 533], [133, 345], [252, 371], [635, 640], [185, 124], [-28, 250], [-387, 279]], [[56639, 89578], [-93, 230], [-8, 910], [-433, 402], [-371, 289]], [[55734, 91409], [167, 156], [309, -312], [362, 29], [298, -143], [265, 262], [137, 433], [431, 200], [356, -235], [-117, -414]], [[99547, 40335], [96, -171], [-46, -308], [-172, -81], [-153, 73], [-27, 260], [107, 203], [126, -74], [69, 98]], [[0, 41087], [57, 27], [-34, -284], [-23, -32], [99822, -145], [-177, -124], [-36, 220], [139, 121], [88, 33], [-99836, 184]], [[33000, 19946], [333, 354], [236, -148], [167, 237], [222, -266], [-83, -207], [-375, -177], [-125, 207], [-236, -266], [-139, 266]], [[34854, 51946], [70, 252], [24, 269], [48, 253], [-107, 349]], [[34889, 53069], [-22, 404], [144, 508]], [[35011, 53981], [95, -65], [204, -140], [294, -499], [46, -242]], [[52655, 75484], [-92, -456], [-126, 120], [-64, 398], [56, 219], [179, 226], [47, -507]], [[51576, 79843], [62, -52], [80, 13]], [[51900, 77752], [-11, -167], [82, -222], [-97, -180], [72, -457], [151, -75], [-32, -256]], [[52065, 76395], [-252, -334], [-548, 160], [-404, -192], [-32, -355]], [[49471, 76235], [144, 354], [53, 1177], [-287, 620], [-205, 299], [-424, 227], [-28, 431], [360, 129], [466, -152], [-88, 669], [263, -254], [646, 461], [84, 484], [243, 119]], [[53081, 48229], [-285, 596], [-184, 488], [-169, 610], [9, 196], [61, 189], [67, 430], [56, 438]], [[52636, 51176], [94, 35], [404, -6], [-2, 711]], [[48278, 82406], [-210, 122], [-172, -9], [57, 317], [-57, 317]], [[47896, 83153], [233, 24], [298, -365], [-149, -406]], [[49165, 85222], [-297, -639], [283, 81], [304, -3], [-72, -481], [-250, -530], [287, -38], [22, -62], [248, -697], [190, -95], [171, -673], [79, -233], [337, -113], [-34, -378], [-142, -173], [111, -305], [-250, -310], [-371, 6], [-473, -163], [-130, 116], [-183, -276], [-257, 67], [-195, -226], [-148, 118], [407, 621], [249, 127], [-2, 1], [-434, 98], [-79, 235], [291, 183], [-152, 319], [52, 387], [413, -54], [1, 0], [40, 343], [-186, 364], [-4, 8], [-337, 104], [-66, 160], [101, 264], [-92, 163], [-149, -279], [-17, 569], [-140, 301], [101, 611], [216, 480], [222, -47], [335, 49]], [[61542, 75120], [42, 252], [-70, 403], [-160, 218], [-154, 68], [-102, 181]], [[61098, 76242], [34, 70], [235, -101], [409, -96], [378, -283], [48, -110], [169, 93], [259, -124], [85, -242], [175, -137]], [[62106, 74858], [-268, 290], [-296, -28]], [[50294, 54083], [-436, -346], [-154, -203], [-250, -171], [-248, 168]], [[50006, 57090], [-20, -184], [116, -305], [-1, -429], [27, -466], [69, -215], [-61, -532], [22, -294], [74, -375], [62, -207]], [[47655, 55121], [-78, 15], [-57, -238], [-78, 3], [-55, 126], [19, 237], [-116, 362], [-73, -67], [-59, -13]], [[47158, 55546], [-77, -34], [3, 217], [-44, 155], [9, 171], [-60, 249], [-78, 211], [-222, 1], [-65, -112], [-76, -13], [-48, -128], [-32, -163], [-148, -260]], [[46320, 55840], [-122, 349], [-108, 232], [-71, 76], [-69, 118], [-32, 261], [-41, 130], [-80, 97]], [[45797, 57103], [123, 288], [84, -11], [73, 99], [61, 1], [44, 78], [-24, 196], [31, 62], [5, 200]], [[46194, 58016], [134, -6], [200, -144], [61, 13], [21, 66], [151, -47], [40, 33]], [[46801, 57931], [16, -216], [44, 1], [73, 78], [46, -19], [77, -150], [119, -48], [76, 128], [90, 79], [67, 83], [55, -15], [62, -130], [33, -163], [114, -248], [-57, -152], [-11, -192], [59, 58], [35, -69], [-15, -176], [85, -170]], [[45321, 58350], [36, 262]], [[45357, 58612], [302, 17], [63, 140], [88, 9], [110, -145], [86, -3], [92, 99], [56, -170], [-120, -133], [-121, 11], [-119, 124], [-103, -136], [-50, -5], [-67, -83], [-253, 13]], [[45797, 57103], [-149, 247], [-117, 39], [-63, 166], [1, 90], [-84, 125], [-18, 127]], [[45367, 57897], [147, 96], [92, -19], [75, 67], [513, -25]], [[52636, 51176], [-52, 90], [96, 663]], [[56583, 71675], [152, -199], [216, 34], [207, -42], [-7, -103], [151, 71], [-35, -175], [-400, -50], [3, 98], [-339, 115], [52, 251]], [[57237, 74699], [-169, 17], [-145, 56], [-336, -154], [192, -332], [-141, -96], [-154, -1], [-147, 305], [-52, -130], [62, -353], [139, -277], [-105, -129], [155, -273], [137, -171], [4, -334], [-257, 157], [82, -302], [-176, -62], [105, -521], [-184, -8], [-228, 257], [-104, 473], [-49, 393], [-108, 272], [-143, 337], [-18, 168]], [[55838, 74710], [182, 53], [106, 129], [150, -12], [46, 103], [53, 20]], [[57254, 75292], [135, -157], [-86, -369], [-66, -67]], [[37010, 99398], [932, 353], [975, -27], [354, 218], [982, 57], [2219, -74], [1737, -469], [-513, -227], [-1062, -26], [-1496, -58], [140, -105], [984, 65], [836, -204], [540, 181], [231, -212], [-305, -344], [707, 220], [1348, 229], [833, -114], [156, -253], [-1132, -420], [-157, -136], [-888, -102], [643, -28], [-324, -431], [-224, -383], [9, -658], [333, -386], [-434, -24], [-457, -187], [513, -313], [65, -502], [-297, -55], [360, -508], [-617, -42], [322, -241], [-91, -208], [-391, -91], [-388, -2], [348, -400], [4, -263], [-549, 244], [-143, -158], [375, -148], [364, -361], [105, -476], [-495, -114], [-214, 228], [-344, 340], [95, -401], [-322, -311], [732, -25], [383, -32], [-745, -515], [-755, -466], [-813, -204], [-306, -2], [-288, -228], [-386, -624], [-597, -414], [-192, -24], [-370, -145], [-399, -138], [-238, -365], [-4, -415], [-141, -388], [-453, -472], [112, -462], [-125, -488], [-142, -577], [-391, -36], [-410, 482], [-556, 3], [-269, 324], [-186, 577], [-481, 735], [-141, 385], [-38, 530], [-384, 546], [100, 435], [-186, 208], [275, 691], [418, 220], [110, 247], [58, 461], [-318, -209], [-151, -88], [-249, -84], [-341, 193], [-19, 401], [109, 314], [258, 9], [567, -157], [-478, 375], [-249, 202], [-276, -83], [-232, 147], [310, 550], [-169, 220], [-220, 409], [-335, 626], [-353, 230], [3, 247], [-745, 346], [-590, 43], [-743, -24], [-677, -44], [-323, 188], [-482, 372], [729, 186], [559, 31], [-1188, 154], [-627, 241], [39, 229], [1051, 285], [1018, 284], [107, 214], [-750, 213], [243, 235], [961, 413], [404, 63], [-115, 265], [658, 156], [854, 93], [853, 5], [303, -184], [737, 325], [663, -221], [390, -46], [577, -192], [-660, 318], [38, 253]], [[24973, 58695], [-142, 103], [-174, 11], [-127, 117], [-149, 244]], [[24381, 59170], [7, 172], [32, 138], [-39, 111], [133, 481], [357, 2], [7, 201], [-45, 36], [-31, 128], [-103, 136], [-103, 198], [125, 1], [1, 333], [259, 1], [257, -7]], [[25297, 59966], [90, -107], [24, 88], [82, -75]], [[25493, 59872], [-127, -225], [-131, -166], [-20, -113], [22, -116], [-58, -150]], [[25179, 59102], [-65, -37], [15, -69], [-52, -66], [-95, -149], [-9, -86]], [[33400, 55523], [183, -217], [171, -385], [8, -304], [105, -14], [149, -289], [109, -205]], [[34125, 54109], [-44, -532], [-169, -154], [15, -139], [-51, -305], [123, -429], [89, -1], [37, -333], [169, -514]], [[33129, 53652], [-188, 448], [75, 163], [-5, 273], [171, 95], [69, 110], [-95, 220], [24, 215], [220, 347]], [[25745, 58251], [-48, 185], [-84, 51]], [[25613, 58487], [19, 237], [-38, 64], [-57, 42], [-122, -70], [-10, 79], [-84, 95], [-60, 118], [-82, 50]], [[25493, 59872], [29, -23], [61, 104], [79, 8], [26, -48], [43, 29], [129, -53], [128, 15], [90, 66], [32, 66], [89, -31], [66, -40], [73, 14], [55, 51], [127, -82], [44, -13], [85, -110], [80, -132], [101, -91], [73, -162]], [[26903, 59440], [-95, 12], [-38, -81], [-97, -77], [-70, 0], [-61, -76], [-56, 27], [-47, 90], [-29, -17], [-36, -141], [-27, 5], [-4, -121], [-97, -163], [-51, -70], [-29, -74], [-82, 120], [-60, -158], [-58, 4], [-65, -14], [6, -290], [-41, -5], [-35, -135], [-86, -25]], [[55230, 77704], [67, -229], [89, -169], [-107, -222]], [[55155, 75778], [-31, -100]], [[55124, 75678], [-261, 218], [-161, 213], [-254, 176], [-233, 434], [56, 45], [-127, 248], [-5, 200], [-179, 93], [-85, -255], [-82, 198], [6, 205], [10, 9]], [[53809, 77462], [194, -20], [51, 100], [94, -97], [109, -11], [-1, 165], [97, 60], [27, 239], [221, 157]], [[54601, 78055], [88, -73], [208, -253], [229, -114], [104, 89]], [[30081, 61241], [-185, 100], [-131, -41], [-169, 43], [-130, -110], [-149, 184], [24, 190], [256, -82], [210, -47], [100, 131], [-127, 256], [2, 226], [-175, 92], [62, 163], [170, -26], [241, -93]], [[54716, 79012], [141, -151], [103, -65], [233, 73], [22, 118], [111, 18], [135, 92], [30, -38], [130, 74], [66, 139], [91, 36], [297, -180], [59, 61]], [[56134, 79189], [155, -161], [19, -159]], [[56308, 78869], [-170, -123], [-131, -401], [-168, -401], [-223, -111]], [[55616, 77833], [-173, 26], [-213, -155]], [[54601, 78055], [-54, 200], [-47, 6]], [[83531, 44530], [-117, -11], [-368, 414], [259, 116], [146, -180], [97, -180], [-17, -159]], [[84713, 45326], [28, -117], [5, -179]], [[84746, 45030], [-181, -441], [-238, -130], [-33, 71], [25, 201], [119, 360], [275, 235]], [[82749, 45797], [100, -158], [172, 48], [69, -251], [-321, -119], [-193, -79], [-149, 5], [95, 340], [153, 5], [74, 209]], [[84139, 45797], [-41, -328], [-417, -168], [-370, 73], [0, 216], [220, 123], [174, -177], [185, 45], [249, 216]], [[80172, 46575], [533, -59], [61, 244], [515, -284], [101, -383], [417, -108], [341, -351], [-317, -225], [-306, 238], [-251, -16], [-288, 44], [-260, 106], [-322, 225], [-204, 59], [-116, -74], [-506, 243], [-48, 254], [-255, 44], [191, 564], [337, -35], [224, -231], [115, -45], [38, -210]], [[87423, 46908], [-143, -402], [-27, 445], [49, 212], [58, 200], [63, -173], [0, -282]], [[85346, 48536], [-104, -196], [-192, 108], [-54, 254], [281, 29], [69, -195]], [[86241, 48752], [101, -452], [-234, 244], [-232, 49], [-157, -39], [-192, 21], [65, 325], [344, 24], [305, -172]], [[89166, 49043], [5, -1925], [4, -1925]], [[89175, 45193], [-247, 485], [-282, 118], [-69, -168], [-352, -18], [118, 481], [175, 164], [-72, 642], [-134, 496], [-538, 500], [-229, 50], [-417, 546], [-82, -287], [-107, -52], [-63, 216], [-1, 257], [-212, 290], [299, 213], [198, -11], [-23, 156], [-407, 1], [-110, 352], [-248, 109], [-117, 293], [374, 143], [142, 192], [446, -242], [44, -220], [78, -955], [287, -354], [232, 627], [319, 356], [247, 1], [238, -206], [206, -212], [298, -113]], [[84788, 51419], [-223, -587], [-209, -113], [-267, 115], [-463, -29], [-243, -85], [-39, -447], [248, -526], [150, 268], [518, 201], [-22, -272], [-121, 86], [-121, -347], [-245, -229], [263, -757], [-50, -203], [249, -682], [-2, -388], [-148, -173], [-109, 207], [134, 484], [-273, -229], [-69, 164], [36, 228], [-200, 346], [21, 576], [-186, -179], [24, -689], [11, -846], [-176, -85], [-119, 173], [79, 544], [-43, 570], [-117, 4], [-86, 405], [115, 387], [40, 469], [139, 891], [58, 243], [237, 439], [217, -174], [350, -82], [319, 25], [275, 429], [48, -132]], [[85746, 51249], [-15, -517], [-143, 58], [-42, -359], [114, -312], [-78, -71], [-112, 374], [-82, 755], [56, 472], [92, 215], [20, -322], [164, -52], [26, -241]], [[80461, 51765], [47, -395], [190, -334], [179, 121], [177, -43], [162, 299], [133, 52], [263, -166], [226, 126], [143, 822], [107, 205], [96, 672], [319, 0], [241, -100]], [[82744, 53024], [-158, -533], [204, -560], [-48, -272], [312, -546], [-329, -70], [-93, -403], [12, -535], [-267, -404], [-7, -589], [-107, -903], [-41, 210], [-316, -266], [-110, 361], [-198, 34], [-139, 189], [-330, -212], [-101, 285], [-182, -32], [-229, 68], [-43, 793], [-138, 164], [-134, 505], [-38, 517], [32, 548], [165, 392]], [[79393, 47122], [-308, -12], [-234, 494], [-356, 482], [-119, 358], [-210, 481], [-138, 443], [-212, 827], [-244, 493], [-81, 508], [-103, 461], [-250, 372], [-145, 506], [-209, 330], [-290, 652], [-24, 300], [178, -24], [430, -114], [246, -577], [215, -401], [153, -246], [263, -635], [283, -9], [233, -405], [161, -495], [211, -270], [-111, -482], [159, -205], [100, -15], [47, -412], [97, -330], [204, -52], [135, -374], [-70, -735], [-11, -914]], [[72530, 68413], [-176, -268], [-108, -553], [269, -224], [262, -289], [362, -332], [381, -76], [160, -301], [215, -56], [334, -138], [231, 10], [32, 234], [-36, 375], [21, 255]], [[77035, 67277], [20, -224], [-97, -108], [23, -364], [-199, 107], [-359, -408], [8, -338], [-153, -496], [-14, -288], [-124, -487], [-217, 135], [-11, -612], [-63, -201], [30, -251], [-137, -140]], [[74730, 63611], [-39, -216], [-189, 7], [-343, -122], [16, -445], [-148, -349], [-400, -398], [-311, -695], [-209, -373], [-276, -387], [-1, -271], [-138, -146], [-251, -212], [-129, -31], [-84, -450], [58, -769], [15, -490], [-118, -561], [-1, -1004], [-144, -29], [-126, -450], [84, -195], [-253, -168], [-93, -401], [-112, -170], [-263, 552], [-128, 827], [-107, 596], [-97, 279], [-148, 568], [-69, 739], [-48, 369], [-253, 811], [-115, 1145], [-83, 756], [1, 716], [-54, 553], [-404, -353], [-196, 70], [-362, 716], [133, 214], [-82, 232], [-326, 501]], [[68937, 64577], [185, 395], [612, -2], [-56, 507], [-156, 300], [-31, 455], [-182, 265], [306, 619], [323, -45], [290, 620], [174, 599], [270, 593], [-4, 421], [236, 342], [-224, 292], [-96, 400], [-99, 517], [137, 255], [421, -144], [310, 88], [268, 496]], [[48278, 82406], [46, -422], [-210, -528], [-493, -349], [-393, 89], [225, 617], [-145, 601], [378, 463], [210, 276]], [[64978, 72558], [244, 114], [197, 338], [186, -17], [122, 110], [197, -55], [308, -299], [221, -65], [318, -523], [207, -21], [24, -498]], [[66909, 68203], [137, -310], [112, -357], [266, -260], [7, -520], [133, -96], [23, -272], [-400, -305], [-105, -687]], [[67082, 65396], [-523, 179], [-303, 136], [-313, 76], [-118, 725], [-133, 105], [-214, -106], [-280, -286], [-339, 196], [-281, 454], [-267, 168], [-186, 561], [-205, 788], [-149, -96], [-177, 196], [-104, -231]], [[63490, 68261], [-153, 311], [-3, 314], [-89, 0], [46, 428], [-143, 449], [-340, 324], [-193, 562], [65, 461], [139, 204], [-21, 345], [-182, 177], [-180, 705]], [[62436, 72541], [-152, 473], [55, 183], [-87, 678], [190, 168]], [[63578, 73220], [88, -436], [263, -123], [193, -296], [395, -102], [434, 156], [27, 139]], [[63490, 68261], [-164, 29]], [[63326, 68290], [-187, 49], [-204, -567]], [[62935, 67772], [-516, 47], [-784, 1188], [-413, 414], [-335, 160]], [[60887, 69581], [-112, 720]], [[60775, 70301], [615, 614], [105, 715], [-26, 431], [152, 146], [142, 369]], [[61763, 72576], [119, 92], [324, -77], [97, -150], [133, 100]], [[45969, 89843], [-64, -382], [314, -403], [-361, -451], [-801, -405], [-240, -107], [-365, 87], [-775, 187], [273, 261], [-605, 289], [492, 114], [-12, 174], [-583, 137], [188, 385], [421, 87], [433, -400], [422, 321], [349, -167], [453, 315], [461, -42]], [[59922, 69905], [-49, -186]], [[59873, 69719], [-100, 82], [-58, -394], [69, -66], [-71, -81], [-12, -156], [131, 80]], [[59832, 69184], [7, -230], [-139, -944]], [[59518, 69025], [80, 194], [-19, 34], [74, 276], [56, 446], [40, 149], [8, 6]], [[59757, 70130], [93, -1], [25, 104], [75, 8]], [[59950, 70241], [4, -242], [-38, -90], [6, -4]], [[54311, 73167], [-100, -465], [41, -183], [-58, -303], [-213, 222], [-141, 64], [-387, 300], [38, 304], [325, -54], [284, 64], [211, 51]], [[52558, 74927], [166, -419], [-39, -782], [-126, 38], [-113, -197], [-105, 156], [-11, 713], [-64, 338], [153, -30], [139, 183]], [[53835, 78058], [-31, -291], [67, -251]], [[53871, 77516], [-221, 86], [-226, -210], [15, -293], [-34, -168], [91, -301], [261, -298], [140, -488], [309, -476], [217, 3], [68, -130], [-78, -118], [249, -214], [204, -178], [238, -308], [29, -111], [-52, -211], [-154, 276], [-242, 97], [-116, -382], [200, -219], [-33, -309], [-116, -35], [-148, -506], [-116, -46], [1, 181], [57, 317], [60, 126], [-108, 342], [-85, 298], [-115, 74], [-82, 255], [-179, 107], [-120, 238], [-206, 38], [-217, 267], [-254, 384], [-189, 340], [-86, 585], [-138, 68], [-226, 195], [-128, -80], [-161, -274], [-115, -43]], [[28453, 61504], [187, -53], [147, -142], [46, -161], [-195, -11], [-84, -99], [-156, 95], [-159, 215], [34, 135], [116, 41], [64, -20]], [[59922, 69905], [309, -234], [544, 630]], [[60887, 69581], [-53, -89], [-556, -296], [277, -591], [-92, -101], [-46, -197], [-212, -82], [-66, -213], [-120, -182], [-310, 94]], [[59709, 67924], [-9, 86]], [[59832, 69184], [41, 173], [0, 362]], [[87399, 70756], [35, -203], [-156, -357], [-114, 189], [-143, -137], [-73, -346], [-181, 168], [2, 281], [154, 352], [158, -68], [114, 248], [204, -127]], [[89159, 72524], [-104, -472], [48, -296], [-145, -416], [-355, -278], [-488, -36], [-396, -675], [-186, 227], [-12, 442], [-483, -130], [-329, -279], [-325, -11], [282, -435], [-186, -1004], [-179, -248], [-135, 229], [69, 533], [-176, 172], [-113, 405], [263, 182], [145, 371], [280, 306], [203, 403], [553, 177], [297, -121], [291, 1050], [185, -282], [408, 591], [158, 229], [174, 723], [-47, 664], [117, 374], [295, 108], [152, -819], [-9, -479], [-256, -595], [4, -610]], [[89974, 76679], [195, -126], [197, 250], [62, -663], [-412, -162], [-244, -587], [-436, 404], [-152, -646], [-308, -9], [-39, 587], [138, 455], [296, 33], [81, 817], [83, 460], [326, -615], [213, -198]], [[69711, 75551], [-159, -109], [-367, -412], [-121, -422], [-104, -4], [-76, 280], [-353, 19], [-57, 484], [-135, 4], [21, 593], [-333, 431], [-476, -46], [-326, -86], [-265, 533], [-227, 223], [-431, 423], [-52, 51], [-715, -349], [11, -2178]], [[65546, 74986], [-142, -29], [-195, 463], [-188, 166], [-315, -123], [-123, -197]], [[64583, 75266], [-15, 144], [68, 246], [-53, 206], [-322, 202], [-125, 530], [-154, 150], [-9, 192], [270, -56], [11, 432], [236, 96], [243, -88], [50, 576], [-50, 365], [-278, -28], [-236, 144], [-321, -260], [-259, -124]], [[63639, 77993], [-142, 96], [29, 304], [-177, 395], [-207, -17], [-235, 401], [160, 448], [-81, 120], [222, 649], [285, -342], [35, 431], [573, 643], [434, 15], [612, -409], [329, -239], [295, 249], [440, 12], [356, -306], [80, 175], [391, -25], [69, 280], [-450, 406], [267, 288], [-52, 161], [266, 153], [-200, 405], [127, 202], [1039, 205], [136, 146], [695, 218], [250, 245], [499, -127], [88, -612], [290, 144], [356, -202], [-23, -322], [267, 33], [696, 558], [-102, -185], [355, -457], [620, -1500], [148, 309], [383, -340], [399, 151], [154, -106], [133, -341], [194, -115], [119, -251], [358, 79], [147, -361]], [[72294, 75601], [-171, 87], [-140, 212], [-412, 62], [-461, 16], [-100, -65], [-396, 248], [-158, -122], [-43, -349], [-457, 204], [-183, -84], [-62, -259]], [[61551, 49585], [-195, -236], [-68, -246], [-104, -44], [-40, -416], [-89, -238], [-54, -393], [-112, -195]], [[60889, 47817], [-399, 590], [-19, 343], [-1007, 1203], [-47, 65]], [[59417, 50018], [-3, 627], [80, 239], [137, 391], [101, 431], [-123, 678], [-32, 296], [-132, 411]], [[59445, 53091], [171, 352], [188, 390]], [[61626, 52895], [-243, -670], [3, -2152], [165, -488]], [[70465, 73876], [-526, -89], [-343, 192], [-301, -46], [26, 340], [303, -98], [101, 182]], [[69725, 74357], [212, -58], [355, 425], [-329, 311], [-198, -147], [-205, 223], [234, 382], [-83, 58]], [[78495, 57780], [-66, 713], [178, 492], [359, 112], [261, -84]], [[79227, 59013], [229, -232], [126, 407], [246, -217]], [[79828, 58971], [64, -394], [-34, -708], [-467, -455], [122, -358], [-292, -43], [-240, -238]], [[78981, 56775], [-233, 87], [-112, 307], [-141, 611]], [[85652, 73393], [240, -697], [68, -383], [3, -681], [-105, -325], [-252, -113], [-222, -245], [-250, -51], [-31, 322], [51, 443], [-122, 615], [206, 99], [-190, 506]], [[85048, 72883], [17, 54], [124, -21], [108, 266], [197, 29], [118, 39], [40, 143]], [[55575, 75742], [52, 132]], [[55627, 75874], [66, 43], [38, 196], [50, 33], [40, -84], [52, -36], [36, -94], [46, -28], [54, -110], [39, 4], [-31, -144], [-33, -71], [9, -44]], [[55993, 75539], [-62, -23], [-164, -91], [-13, -121], [-35, 5]], [[63326, 68290], [58, -261], [-25, -135], [89, -445]], [[63448, 67449], [-196, -16], [-69, 282], [-248, 57]], [[79227, 59013], [90, 266], [12, 500], [-224, 515], [-18, 583], [-211, 480], [-210, 40], [-56, -205], [-163, -17], [-83, 104], [-293, -353], [-6, 530], [68, 623], [-188, 27], [-16, 355], [-120, 182]], [[77809, 62643], [59, 218], [237, 384]], [[78380, 63852], [162, -466], [125, -537], [342, -5], [108, -515], [-178, -155], [-80, -212], [333, -353], [231, -699], [175, -520], [210, -411], [70, -418], [-50, -590]], [[59757, 70130], [99, 482], [138, 416], [5, 21]], [[59999, 71049], [125, -31], [45, -231], [-151, -223], [-68, -323]], [[47857, 53158], [-73, -5], [-286, 282], [-252, 449], [-237, 324], [-187, 381]], [[46822, 54589], [66, 189], [15, 172], [126, 320], [129, 276]], [[54125, 64088], [-197, -220], [-156, 324], [-439, 255]], [[52633, 68486], [136, 137], [24, 250], [-30, 244], [191, 228], [86, 189], [135, 170], [16, 454]], [[53191, 70158], [326, -204], [117, 51], [232, -98], [368, -264], [130, -526], [250, -114], [391, -248], [296, -293], [136, 153], [133, 272], [-65, 452], [87, 288], [200, 277], [192, 80], [375, -121], [95, -264], [104, -2], [88, -101], [276, -70], [68, -195]], [[56944, 63578], [0, -1180], [-320, -2], [-3, -248]], [[56621, 62148], [-1108, 1131], [-1108, 1132], [-280, -323]], [[72718, 55024], [-42, -615], [-116, -168], [-242, -135], [-132, 470], [-49, 849], [126, 959], [192, -328], [129, -416], [134, -616]], [[58049, 33472], [96, -178], [-85, -288], [-47, -192], [-155, -93], [-51, -188], [-99, -59], [-209, 454], [148, 374], [151, 232], [130, 120], [121, -182]], [[56314, 82678], [-23, 150], [30, 162], [-123, 94], [-291, 103]], [[55907, 83187], [-59, 497]], [[55848, 83684], [318, 181], [466, -38], [273, 59], [39, -123], [148, -38], [267, -287]], [[56523, 82432], [-67, 182], [-142, 64]], [[55848, 83684], [10, 445], [136, 371], [262, 202], [221, -442], [223, 12], [53, 453]], [[57579, 84537], [134, -136], [24, -287], [89, -348]], [[47592, 66920], [-42, 0], [7, -317], [-172, -19], [-90, -134], [-126, 0], [-100, 76], [-234, -63], [-91, -460], [-86, -44], [-131, -745], [-386, -637], [-92, -816], [-114, -265], [-33, -213], [-625, -48], [-5, 1]], [[45272, 63236], [13, 274], [106, 161], [91, 308], [-18, 200], [96, 417], [155, 376], [93, 95], [74, 344], [6, 315], [100, 365], [185, 216], [177, 603], [5, 8], [139, 227], [259, 65], [218, 404], [140, 158], [232, 493], [-70, 735], [106, 508], [37, 312], [179, 399], [278, 270], [206, 244], [186, 612], [87, 362], [205, -2], [167, -251], [264, 41], [288, -131], [121, -6]], [[57394, 79070], [66, 87], [185, 58], [204, -184], [115, -22], [125, -159], [-20, -200], [101, -97], [40, -247], [97, -150], [-19, -88], [52, -60], [-74, -44], [-164, 18], [-27, 81], [-58, -47], [20, -106], [-76, -188], [-49, -203], [-70, -64]], [[57842, 77455], [-50, 270], [30, 252], [-9, 259], [-160, 352], [-89, 249], [-86, 175], [-84, 58]], [[63761, 43212], [74, -251], [69, -390], [45, -711], [72, -276], [-28, -284], [-49, -174], [-94, 347], [-53, -175], [53, -438], [-24, -250], [-77, -137], [-18, -500], [-109, -689], [-137, -814], [-172, -1120], [-106, -821], [-125, -685], [-226, -140], [-243, -250], [-160, 151], [-220, 211], [-77, 312], [-18, 524], [-98, 471], [-26, 425], [50, 426], [128, 102], [1, 197], [133, 447], [25, 377], [-65, 280], [-52, 372], [-23, 544], [97, 331], [38, 375], [138, 22], [155, 121], [103, 107], [122, 7], [158, 337], [229, 364], [83, 297], [-38, 253], [118, -71], [153, 410], [6, 356], [92, 264], [96, -254]], [[23016, 65864], [-107, -518], [-49, -426], [-20, -791], [-27, -289], [48, -322], [86, -288], [56, -458], [184, -440], [65, -337], [109, -291], [295, -157], [114, -247], [244, 165], [212, 60], [208, 106], [175, 101], [176, 241], [67, 345], [22, 496], [48, 173], [188, 155], [294, 137], [246, -21], [169, 50], [66, -125], [-9, -285], [-149, -351], [-66, -360], [51, -103], [-42, -255], [-69, -461], [-71, 152], [-58, -10]], [[24381, 59170], [-314, 636], [-144, 191], [-226, 155], [-156, -43], [-223, -223], [-140, -58], [-196, 156], [-208, 112], [-260, 271], [-208, 83], [-314, 275], [-233, 282], [-70, 158], [-155, 35], [-284, 187], [-116, 270], [-299, 335], [-139, 373], [-66, 288], [93, 57], [-29, 169], [64, 153], [1, 204], [-93, 266], [-25, 235], [-94, 298], [-244, 587], [-280, 462], [-135, 368], [-238, 241], [-51, 145], [42, 365], [-142, 138], [-164, 287], [-69, 412], [-149, 48], [-162, 311], [-130, 288], [-12, 184], [-149, 446], [-99, 452], [5, 227], [-201, 234], [-93, -25], [-159, 163], [-44, -240], [46, -284], [27, -444], [95, -243], [206, -407], [46, -139], [42, -42], [37, -203], [49, 8], [56, -381], [85, -150], [59, -210], [174, -300], [92, -550], [83, -259], [77, -277], [15, -311], [134, -20], [112, -268], [100, -264], [-6, -106], [-117, -217], [-49, 3], [-74, 359], [-181, 337], [-201, 286], [-142, 150], [9, 432], [-42, 320], [-132, 183], [-191, 264], [-37, -76], [-70, 154], [-171, 143], [-164, 343], [20, 44], [115, -33], [103, 221], [10, 266], [-214, 422], [-163, 163], [-102, 369], [-103, 388], [-129, 472], [-113, 531]], [[17464, 69802], [316, 46], [353, 64], [-26, -116], [419, -287], [634, -416], [552, 4], [221, 0], [0, 244], [481, 0], [102, -210], [142, -186], [165, -260], [92, -309], [69, -325], [144, -178], [230, -177], [175, 467], [227, 11], [196, -236], [139, -404], [96, -346], [164, -337], [61, -414], [78, -277], [217, -184], [197, -130], [108, 18]], [[55993, 75539], [95, 35], [128, 9]], [[46619, 59216], [93, 107], [47, 348], [88, 14], [194, -165], [157, 117], [107, -39], [42, 131], [1114, 9], [62, 414], [-48, 73], [-134, 2550], [-134, 2550], [425, 10]], [[51185, 61897], [1, -1361], [-152, -394], [-24, -364], [-247, -94], [-379, -51], [-102, -210], [-178, -23]], [[46801, 57931], [13, 184], [-24, 229], [-104, 166], [-54, 338], [-13, 368]], [[77375, 56448], [-27, 439], [86, 452], [-94, 350], [23, 644], [-113, 306], [-90, 707], [-50, 746], [-121, 490], [-183, -297], [-315, -421], [-156, 53], [-172, 138], [96, 732], [-58, 554], [-218, 681], [34, 213], [-163, 76], [-197, 481]], [[77809, 62643], [-159, -137], [-162, -256], [-196, -26], [-127, -639], [-117, -107], [134, -519], [177, -431], [113, -390], [-101, -514], [-96, -109], [66, -296], [185, -470], [32, -330], [-4, -274], [108, -539], [-152, -551], [-135, -607]], [[55380, 75322], [-58, 46], [-78, 192], [-120, 118]], [[55338, 76294], [74, -101], [40, -82], [91, -63], [106, -123], [-22, -51]], [[74375, 79706], [292, 102], [530, 509], [423, 278], [242, -182], [289, -8], [186, -276], [277, -22], [402, -148], [270, 411], [-113, 348], [288, 612], [311, -244], [252, -69], [327, -152], [53, -443], [394, -248], [263, 109], [351, 78], [279, -78], [272, -284], [168, -302], [258, 6], [350, -96], [255, 146], [366, 98], [407, 416], [166, -63], [146, -198], [331, 49]], [[59599, 43773], [209, 48], [334, -166], [73, 74], [193, 16], [99, 177], [167, -10], [303, 230], [221, 342]], [[61198, 44484], [45, -265], [-11, -588], [34, -519], [11, -923], [49, -290], [-83, -422], [-108, -410], [-177, -366], [-254, -225], [-313, -287], [-313, -634], [-107, -108], [-194, -420], [-115, -136], [-23, -421], [132, -448], [54, -346], [4, -177], [49, 29], [-8, -579], [-45, -275], [65, -101], [-41, -245], [-116, -211], [-229, -199], [-334, -320], [-122, -219], [24, -248], [71, -40], [-24, -311]], [[59119, 34780], [-211, 5]], [[58908, 34785], [-24, 261], [-41, 265]], [[58843, 35311], [-23, 212], [49, 659], [-72, 419], [-133, 832]], [[58664, 37433], [292, 671], [74, 426], [42, 53], [31, 348], [-45, 175], [12, 442], [54, 409], [0, 748], [-145, 190], [-132, 43], [-60, 146], [-128, 125], [-232, -12], [-18, 220]], [[58409, 41417], [-26, 421], [843, 487]], [[59226, 42325], [159, -284], [77, 54], [110, -149], [16, -237], [-59, -274], [21, -417], [181, -365], [85, 410], [120, 124], [-24, 760], [-116, 427], [-100, 191], [-97, -9], [-77, 768], [77, 449]], [[46619, 59216], [-184, 405], [-168, 435], [-184, 157], [-133, 173], [-155, -6], [-135, -129], [-138, 51], [-96, -189]], [[45426, 60113], [-24, 318], [78, 291], [34, 557], [-30, 583], [-34, 294], [28, 295], [-72, 281], [-146, 255]], [[45260, 62987], [60, 197], [1088, -4], [-53, 853], [68, 304], [261, 53], [-9, 1512], [911, -31], [1, 895]], [[59226, 42325], [-147, 153], [85, 549], [87, 205], [-53, 490], [56, 479], [47, 160], [-71, 501], [-131, 264]], [[59099, 45126], [273, -110], [55, -164], [95, -275], [77, -804]], [[78372, 54256], [64, -56], [164, -356], [116, -396], [16, -398], [-29, -269], [27, -203], [20, -349], [98, -163], [109, -523], [-5, -199], [-197, -40], [-263, 438], [-329, 469], [-32, 301], [-161, 395], [-38, 489], [-100, 322], [30, 431], [-61, 250]], [[77801, 54399], [48, 105], [227, -258], [22, -304], [183, 71], [91, 243]], [[80461, 51765], [204, -202], [214, 110], [56, 500], [119, 112], [333, 128], [199, 467], [137, 374]], [[82069, 53798], [214, 411], [140, 462], [112, 2], [143, -299], [13, -257], [183, -165], [231, -177], [-20, -232], [-186, -29], [50, -289], [-205, -201]], [[54540, 33696], [-207, 446], [-108, 432], [-62, 575], [-68, 428], [-93, 910], [-7, 707], [-35, 322], [-108, 243], [-144, 489], [-146, 708], [-60, 371], [-226, 577], [-17, 453]], [[56448, 40227], [228, 134], [180, -34], [109, -133], [2, -49]], [[55526, 35946], [0, -2182], [-248, -302], [-149, -43], [-175, 112], [-125, 43], [-47, 252], [-109, 162], [-133, -292]], [[96049, 38125], [228, -366], [144, -272], [-105, -142], [-153, 160], [-199, 266], [-179, 313], [-184, 416], [-38, 201], [119, -9], [156, -201], [122, -200], [89, -166]], [[54125, 64088], [68, -919], [104, -153], [4, -188], [116, -203], [-60, -254], [-107, -1199], [-15, -769], [-354, -557], [-120, -778], [115, -219], [0, -380], [178, -13], [-28, -279]], [[53939, 57955], [-52, -13], [-188, 647], [-65, 24], [-217, -331], [-215, 173], [-150, 34], [-80, -83], [-163, 18], [-164, -252], [-141, -14], [-337, 305], [-131, -145], [-142, 10], [-104, 223], [-279, 221], [-298, -70], [-72, -128], [-39, -340], [-80, -238], [-19, -527]], [[52361, 53399], [-289, -213], [-105, 31], [-107, -132], [-222, 13], [-149, 370], [-91, 427], [-197, 389], [-209, -7], [-245, 1]], [[26191, 57131], [-96, 186], [-130, 238], [-61, 200], [-117, 185], [-140, 267], [31, 91], [46, -88], [21, 41]], [[26903, 59440], [-24, -57], [-14, -132], [29, -216], [-64, -202], [-30, -237], [-9, -261], [15, -152], [7, -266], [-43, -58], [-26, -253], [19, -156], [-56, -151], [12, -159], [43, -97]], [[50920, 80916], [143, 162], [244, 869], [380, 248], [231, -17]], [[58639, 91676], [-473, -237], [-224, -54]], [[55734, 91409], [-172, -24], [-41, -389], [-523, 95], [-74, -329], [-267, 2], [-183, -421], [-278, -655], [-431, -831], [101, -202], [-97, -234], [-275, 10], [-180, -554], [17, -784], [177, -300], [-92, -694], [-231, -405], [-122, -341]], [[53063, 85353], [-187, 363], [-548, -684], [-371, -138], [-384, 301], [-99, 635], [-88, 1363], [256, 381], [733, 496], [549, 609], [508, 824], [668, 1141], [465, 444], [763, 741], [610, 259], [457, -31], [423, 489], [506, -26], [499, 118], [869, -433], [-358, -158], [305, -371]], [[56867, 96577], [-620, -241], [-490, 137], [191, 152], [-167, 189], [575, 119], [110, -222], [401, -134]], [[55069, 97669], [915, -440], [-699, -233], [-155, -435], [-243, -111], [-132, -490], [-335, -23], [-598, 361], [252, 210], [-416, 170], [-541, 499], [-216, 463], [757, 212], [152, -207], [396, 8], [105, 202], [408, 20], [350, -206]], [[57068, 98086], [545, -207], [-412, -318], [-806, -70], [-819, 98], [-50, 163], [-398, 11], [-304, 271], [858, 165], [403, -142], [281, 177], [702, -148]], [[98060, 26404], [63, -244], [198, 239], [80, -249], [0, -249], [-103, -274], [-182, -435], [-142, -238], [103, -284], [-214, -7], [-238, -223], [-75, -387], [-157, -597], [-219, -264], [-138, -169], [-256, 13], [-180, 194], [-302, 42], [-46, 217], [149, 438], [349, 583], [179, 111], [200, 225], [238, 310], [167, 306], [123, 441], [106, 149], [41, 330], [195, 273], [61, -251]], [[98502, 29218], [202, -622], [5, 403], [126, -161], [41, -447], [224, -192], [188, -48], [158, 226], [141, -69], [-67, -524], [-85, -345], [-212, 12], [-74, -179], [26, -254], [-41, -110], [-105, -319], [-138, -404], [-214, -236], [-48, 155], [-116, 85], [160, 486], [-91, 326], [-299, 236], [8, 214], [201, 206], [47, 455], [-13, 382], [-113, 396], [8, 104], [-133, 244], [-218, 523], [-117, 418], [104, 46], [151, -328], [216, -153], [78, -526]], [[64752, 60417], [-91, 413], [-217, 975]], [[64444, 61805], [833, 591], [185, 1182], [-127, 418]], [[65665, 65306], [125, -404], [155, -214], [203, -78], [165, -107], [125, -339], [75, -196], [100, -75], [-1, -132], [-101, -352], [-44, -166], [-117, -189], [-104, -404], [-126, 31], [-58, -141], [-44, -300], [34, -395], [-26, -72], [-128, 2], [-174, -221], [-27, -288], [-63, -125], [-173, 5], [-109, -149], [1, -238], [-134, -165], [-153, 56], [-186, -199], [-128, -34]], [[65575, 65974], [80, 201], [35, -51], [-26, -244], [-37, -108]], [[68937, 64577], [-203, 150], [-83, 424], [-215, 450], [-512, -111], [-451, -11], [-391, -83]], [[28366, 54848], [-93, 170], [-59, 319], [68, 158], [-70, 40], [-52, 196], [-138, 164], [-122, -38], [-56, -205], [-112, -149], [-61, -20], [-27, -123], [132, -321], [-75, -76], [-40, -87], [-130, -30], [-48, 353], [-36, -101], [-92, 35], [-56, 238], [-114, 39], [-72, 69], [-119, -1], [-8, -128], [-32, 89]], [[27070, 56232], [100, -212], [-6, -126], [111, -26], [26, 48], [77, -145], [136, 42], [119, 150], [168, 119], [95, 176], [153, -34], [-10, -58], [155, -21], [124, -102], [90, -177], [105, -164]], [[30452, 39739], [-279, 340], [-24, 242], [-551, 593], [-498, 646], [-214, 365], [-115, 488], [46, 170], [-236, 775], [-274, 1090], [-262, 1177], [-114, 269], [-87, 435], [-216, 386], [-198, 239], [90, 264], [-134, 563], [86, 414], [221, 373]], [[85104, 55551], [28, -392], [16, -332], [-94, -540], [-102, 602], [-130, -300], [89, -435], [-79, -277], [-327, 343], [-78, 428], [84, 280], [-176, 280], [-87, -245], [-131, 23], [-205, -330], [-46, 173], [109, 498], [175, 166], [151, 223], [98, -268], [212, 162], [45, 264], [196, 15], [-16, 457], [225, -280], [23, -297], [20, -218]], [[84439, 56653], [-100, -195], [-87, -373], [-87, -175], [-171, 409], [57, 158], [70, 165], [30, 367], [153, 35], [-44, -398], [205, 570], [-26, -563]], [[82917, 56084], [-369, -561], [136, 414], [200, 364], [167, 409], [146, 587], [49, -482], [-183, -325], [-146, -406]], [[83856, 57606], [166, -183], [177, 1], [-5, -247], [-129, -251], [-176, -178], [-10, 275], [20, 301], [-43, 282]], [[84861, 57766], [78, -660], [-214, 157], [5, -199], [68, -364], [-132, -133], [-11, 416], [-84, 31], [-43, 357], [163, -47], [-4, 224], [-169, 451], [266, -13], [77, -220]], [[83757, 58301], [-74, -510], [-119, 295], [-142, 450], [238, -22], [97, -213]], [[83700, 61512], [171, -168], [85, 153], [26, -150], [-46, -245], [95, -423], [-73, -491], [-164, -196], [-43, -476], [62, -471], [147, -65], [123, 70], [347, -328], [-27, -321], [91, -142], [-29, -272], [-216, 290], [-103, 310], [-71, -217], [-177, 354], [-253, -87], [-138, 130], [14, 244], [87, 151], [-83, 136], [-36, -213], [-137, 340], [-41, 257], [-11, 566], [112, -195], [29, 925], [90, 535], [169, -1]], [[93299, 46550], [-78, -59], [-120, 227], [-122, 375], [-59, 450], [38, 57], [30, -175], [84, -134], [135, -375], [131, -200], [-39, -166]], [[92217, 47343], [-146, -48], [-44, -166], [-152, -144], [-142, -138], [-148, 1], [-228, 171], [-158, 165], [23, 183], [249, -86], [152, 46], [42, 283], [40, 15], [27, -314], [158, 45], [78, 202], [155, 211], [-30, 348], [166, 11], [56, -97], [-5, -327], [-93, -361]], [[89166, 49043], [482, -407], [513, -338], [192, -302], [154, -297], [43, -349], [462, -365], [68, -313], [-256, -64], [62, -393], [248, -388], [180, -627], [159, 20], [-11, -262], [215, -100], [-84, -111], [295, -249], [-30, -171], [-184, -41], [-69, 153], [-238, 66], [-281, 89], [-216, 377], [-158, 325], [-144, 517], [-362, 259], [-235, -169], [-170, -195], [35, -436], [-218, -203], [-155, 99], [-288, 25]], [[92538, 47921], [-87, -157], [-52, 348], [-65, 229], [-126, 193], [-158, 252], [-200, 174], [77, 143], [150, -166], [94, -130], [117, -142], [111, -248], [106, -189], [33, -307]], [[53922, 82340], [189, 174], [434, 273], [350, 200], [277, -100], [21, -144], [268, -7]], [[55461, 82736], [342, -67], [511, 9]], [[56535, 81053], [139, -515], [-29, -166], [-138, -69], [-252, -491], [71, -266], [-60, 35]], [[56266, 79581], [-264, 227], [-200, -84], [-131, 61], [-165, -127], [-140, 210], [-114, -81], [-16, 36]], [[31588, 61519], [142, -52], [50, -118], [-71, -149], [-209, 4], [-163, -21], [-16, 253], [40, 86], [227, -3]], [[86288, 75628], [39, -104]], [[86327, 75524], [-106, 36], [-120, -200], [-83, -202], [10, -424], [-143, -130], [-50, -105], [-104, -174], [-185, -97], [-121, -159], [-9, -256], [-32, -65], [111, -96], [157, -259]], [[85048, 72883], [-135, 112], [-34, -111], [-81, -49], [-10, 112], [-72, 54], [-75, 94], [76, 260], [66, 69], [-25, 108], [71, 319], [-18, 96], [-163, 65], [-131, 158]], [[47929, 72498], [-112, -153], [-146, 83], [-143, -65], [42, 462], [-26, 363], [-124, 55], [-67, 224], [22, 386], [111, 215], [20, 239], [58, 355], [-6, 250], [-56, 212], [-12, 200]], [[64113, 65205], [-18, 430], [75, 310], [76, 64], [84, -185], [5, -346], [-61, -348]], [[64274, 65130], [-77, -42], [-84, 117]], [[56308, 78869], [120, 127], [172, -65], [178, -3], [129, -144], [95, 91], [205, 56], [69, 139], [118, 0]], [[57842, 77455], [124, -109], [131, 95], [126, -101]], [[58223, 77340], [6, -152], [-135, -128], [-84, 56], [-78, -713]], [[56293, 76715], [-51, 103], [65, 99], [-69, 74], [-87, -133], [-162, 172], [-22, 244], [-169, 139], [-31, 188], [-151, 232]], [[89901, 80562], [280, -1046], [-411, 195], [-171, -854], [271, -605], [-8, -413], [-211, 356], [-182, -457], [-51, 496], [31, 575], [-32, 638], [64, 446], [13, 790], [-163, 581], [24, 808], [257, 271], [-110, 274], [123, 83], [73, -391], [96, -569], [-7, -581], [114, -597]], [[55461, 82736], [63, 260], [383, 191]], [[99999, 92429], [-305, -30], [-49, 187], [-99645, 247], [36, 24], [235, -1], [402, -169], [-24, -81], [-286, -141], [-363, -36], [99999, 0]], [[89889, 93835], [-421, -4], [-569, 66], [-49, 31], [263, 234], [348, 54], [394, -226], [34, -155]], [[91869, 94941], [-321, -234], [-444, 53], [-516, 233], [66, 192], [518, -89], [697, -155]], [[90301, 95224], [-219, -439], [-1023, 16], [-461, -139], [-550, 384], [149, 406], [366, 111], [734, -26], [1004, -313]], [[65981, 92363], [-164, -52], [-907, 77], [-74, 262], [-503, 158], [-40, 320], [284, 126], [-10, 323], [551, 503], [-255, 73], [665, 518], [-75, 268], [621, 312], [917, 380], [925, 110], [475, 220], [541, 76], [193, -233], [-187, -184], [-984, -293], [-848, -282], [-863, -562], [-414, -577], [-435, -568], [56, -491], [531, -484]], [[63639, 77993], [-127, -350], [-269, -97], [-276, -610], [252, -561], [-27, -398], [303, -696]], [[61098, 76242], [-354, 499], [-317, 223], [-240, 347], [202, 95], [231, 494], [-156, 234], [410, 241], [-8, 129], [-249, -95]], [[60617, 78409], [9, 262], [143, 165], [269, 43], [44, 197], [-62, 326], [113, 310], [-3, 173], [-410, 192], [-162, -6], [-172, 277], [-213, -94], [-352, 208], [6, 116], [-99, 256], [-222, 29], [-23, 183], [70, 120], [-178, 334], [-288, -57], [-84, 30], [-70, -134], [-104, 23]], [[57772, 85719], [316, 327], [-291, 280]], [[58639, 91676], [286, 206], [456, -358], [761, -140], [1050, -668], [213, -281], [18, -393], [-308, -311], [-454, -157], [-1240, 449], [-204, -75], [453, -433], [18, -274], [18, -604], [358, -180], [217, -153], [36, 286], [-168, 254], [177, 224], [672, -368], [233, 144], [-186, 433], [647, 578], [256, -34], [260, -206], [161, 406], [-231, 352], [136, 353], [-204, 367], [777, -190], [158, -331], [-351, -73], [1, -328], [219, -203], [429, 128], [68, 377], [580, 282], [970, 507], [209, -29], [-273, -359], [344, -61], [199, 202], [521, 16], [412, 245], [317, -356], [315, 391], [-291, 343], [145, 195], [820, -179], [385, -185], [1006, -675], [186, 309], [-282, 313], [-8, 125], [-335, 58], [92, 280], [-149, 461], [-8, 189], [512, 535], [183, 537], [206, 116], [736, -156], [57, -328], [-263, -479], [173, -189], [89, -413], [-63, -809], [307, -362], [-120, -395], [-544, -839], [318, -87], [110, 213], [306, 151], [74, 293], [240, 281], [-162, 336], [130, 390], [-304, 49], [-67, 328], [222, 593], [-361, 482], [497, 398], [-64, 421], [139, 13], [145, -328], [-109, -570], [297, -108], [-127, 426], [465, 233], [577, 31], [513, -337], [-247, 492], [-28, 630], [483, 119], [669, -26], [602, 77], [-226, 309], [321, 388], [319, 16], [540, 293], [734, 79], [93, 162], [729, 55], [227, -133], [624, 314], [510, -10], [77, 255], [265, 252], [656, 242], [476, -191], [-378, -146], [629, -90], [75, -292], [254, 143], [812, -7], [626, -289], [223, -221], [-69, -307], [-307, -175], [-730, -328], [-209, -175], [345, -83], [410, -149], [251, 112], [141, -379], [122, 153], [444, 93], [892, -97], [67, -276], [1162, -88], [15, 451], [590, -104], [443, 4], [449, -312], [128, -378], [-165, -247], [349, -465], [437, -240], [268, 620], [446, -266], [473, 159], [538, -182], [204, 166], [455, -83], [-201, 549], [367, 256], [2509, -384], [236, -351], [727, -451], [1122, 112], [553, -98], [231, -244], [-33, -432], [342, -168], [372, 121], [492, 15], [525, -116], [526, 66], [484, -526], [344, 189], [-224, 378], [123, 262], [886, -165], [578, 36], [799, -282], [-99610, -258], [681, -451], [728, -588], [-24, -367], [187, -147], [-64, 429], [754, -88], [544, -553], [-276, -257], [-455, -61], [-7, -578], [-111, -122], [-260, 17], [-212, 206], [-369, 172], [-62, 257], [-283, 96], [-315, -76], [-151, 207], [60, 219], [-333, -140], [126, -278], [-158, -251], [99997, -3], [-357, -260], [-360, 44], [250, -315], [166, -487], [128, -159], [32, -244], [-71, -157], [-518, 129], [-777, -445], [-247, -69], [-425, -415], [-403, -362], [-102, -269], [-397, 409], [-724, -464], [-126, 219], [-268, -253], [-371, 81], [-90, -388], [-333, -572], [10, -239], [316, -132], [-37, -860], [-258, -22], [-119, -494], [116, -255], [-486, -302], [-96, -674], [-415, -144], [-83, -600], [-400, -551], [-103, 407], [-119, 862], [-155, 1313], [134, 819], [234, 353], [14, 276], [432, 132], [496, 744], [479, 608], [499, 471], [223, 833], [-337, -50], [-167, -487], [-705, -649], [-227, 727], [-717, -201], [-696, -990], [230, -362], [-620, -154], [-430, -61], [20, 427], [-431, 90], [-344, -291], [-850, 102], [-914, -175], [-899, -1153], [-1065, -1394], [438, -74], [136, -370], [270, -132], [178, 295], [305, -38], [401, -650], [9, -503], [-217, -590], [-23, -705], [-126, -945], [-418, -855], [-94, -409], [-377, -688], [-374, -682], [-179, -349], [-370, -346], [-175, -8], [-175, 287], [-373, -432], [-43, -197]], [[79187, 96845], [-1566, -228], [507, 776], [229, 66], [208, -38], [704, -336], [-82, -240]], [[64204, 98169], [-373, -78], [-250, -45], [-39, -97], [-324, -98], [-301, 140], [158, 185], [-618, 18], [542, 107], [422, 8], [57, -160], [159, 142], [262, 97], [412, -129], [-107, -90]], [[77760, 97184], [-606, -73], [-773, 170], [-462, 226], [-213, 423], [-379, 117], [722, 404], [600, 133], [540, -297], [640, -572], [-69, -531]], [[58449, 49909], [110, -333], [-16, -348], [-80, -74]], [[58216, 49787], [67, -60], [166, 182]], [[45260, 62987], [12, 249]], [[61883, 60238], [-37, 252], [-83, 178], [-22, 236], [-143, 212], [-148, 495], [-79, 482], [-192, 406], [-124, 97], [-184, 563], [-32, 411], [12, 350], [-159, 655], [-130, 231], [-150, 122], [-92, 339], [15, 133], [-77, 306], [-81, 132], [-108, 440], [-170, 476], [-141, 406], [-139, -3], [44, 325], [12, 206], [34, 236]], [[63448, 67449], [109, -510], [137, -135], [47, -207], [190, -249], [16, -243], [-27, -197], [35, -199], [80, -165], [37, -194], [41, -145]], [[64274, 65130], [53, -226]], [[64444, 61805], [-801, -226], [-259, -266], [-199, -620], [-130, -99], [-70, 197], [-106, -30], [-269, 60], [-50, 59], [-321, -14], [-75, -53], [-114, 153], [-74, -290], [28, -249], [-121, -189]], [[59434, 56171], [-39, 12], [5, 294], [-33, 203], [-143, 233], [-34, 426], [34, 436], [-129, 41], [-19, -132], [-167, -30], [67, -173], [23, -355], [-152, -324], [-138, -426], [-144, -61], [-233, 345], [-105, -122], [-29, -172], [-143, -112], [-9, -122], [-277, 0], [-38, 122], [-200, 20], [-100, -101], [-77, 51], [-143, 344], [-48, 163], [-200, -81], [-76, -274], [-72, -528], [-95, -111], [-85, -65]], [[56635, 55672], [-23, 28]], [[56351, 57163], [3, 143], [-102, 174], [-3, 343], [-58, 228], [-98, -34], [28, 217], [72, 246], [-32, 245], [92, 181], [-58, 138], [73, 365], [127, 435], [240, -41], [-14, 2345]], [[60240, 63578], [90, -580], [-61, -107], [40, -608], [102, -706], [106, -145], [152, -219]], [[59433, 56242], [1, -71]], [[59434, 56171], [3, -460]], [[59445, 53091], [-171, -272], [-195, 1], [-224, -138], [-176, 132], [-115, -161]], [[56824, 55442], [-189, 230]], [[45357, 58612], [-115, 460], [-138, 210], [122, 112], [134, 415], [66, 304]], [[45367, 57897], [-46, 453]], [[95032, 44386], [78, -203], [-194, 4], [-106, 363], [166, -142], [56, -22]], [[94680, 44747], [-108, -14], [-170, 60], [-58, 91], [17, 235], [183, -93], [91, -124], [45, -155]], [[94910, 44908], [-42, -109], [-206, 512], [-57, 353], [94, 0], [100, -473], [111, -283]], [[94409, 45654], [12, -119], [-218, 251], [-152, 212], [-104, 197], [41, 60], [128, -142], [228, -272], [65, -187]], [[93760, 46238], [-56, -33], [-121, 134], [-114, 243], [14, 99], [166, -250], [111, -193]], [[46822, 54589], [-75, 44], [-200, 238], [-144, 316], [-49, 216], [-34, 437]], [[25613, 58487], [-31, -139], [-161, 9], [-100, 57], [-115, 117], [-154, 37], [-79, 127]], [[61984, 57352], [91, -109], [54, -245], [125, -247], [138, -2], [262, 151], [302, 70], [245, 184], [138, 39], [99, 108], [158, 20]], [[63596, 57321], [-2, -9], [-1, -244], [0, -596], [0, -308], [-125, -363], [-194, -493]], [[63596, 57321], [89, 12], [128, 88], [147, 59], [132, 202], [105, 2], [6, -163], [-25, -344], [1, -310], [-59, -214], [-78, -639], [-134, -659], [-172, -755], [-238, -866], [-237, -661], [-327, -806], [-278, -479], [-415, -586], [-259, -450], [-304, -715], [-64, -312], [-63, -140]], [[34125, 54109], [333, -119], [30, 107], [225, 43], [298, -159]], [[34889, 53069], [109, -351], [-49, -254], [-24, -270], [-71, -248]], [[56266, 79581], [-77, -154], [-55, -238]], [[53809, 77462], [62, 54]], [[56639, 89578], [-478, -167], [-269, -413], [43, -361], [-441, -475], [-537, -509], [-202, -832], [198, -416], [265, -328], [-255, -666], [-289, -138], [-106, -992], [-157, -554], [-337, 57], [-158, -468], [-321, -27], [-89, 558], [-232, 671], [-211, 835]], [[58908, 34785], [-56, -263], [-163, -63], [-166, 320], [-2, 204], [76, 222], [26, 172], [80, 42], [140, -108]], [[59999, 71049], [-26, 452], [68, 243]], [[60041, 71744], [74, 129], [75, 130], [15, 329], [91, -115], [306, 165], [147, -112], [229, 2], [320, 222], [149, -10], [316, 92]], [[50518, 54209], [-224, -126]], [[78495, 57780], [-249, 271], [-238, -11], [41, 464], [-245, -3], [-22, -650], [-150, -863], [-90, -522], [19, -428], [181, -18], [113, -539], [50, -512], [155, -338], [168, -69], [144, -306]], [[77801, 54399], [-110, 227], [-47, 292], [-148, 334], [-135, 280], [-45, -347], [-53, 328], [30, 369], [82, 566]], [[68841, 72526], [156, 598], [-60, 440], [-204, 140], [72, 261], [232, -28], [132, 326], [89, 380], [371, 137], [-58, -274], [40, -164], [114, 15]], [[64978, 72558], [-52, 417], [40, 618], [-216, 200], [71, 405], [-184, 34], [61, 498], [262, -145], [244, 189], [-202, 355], [-80, 338], [-224, -151], [-28, -433], [-87, 383]], [[65546, 74986], [313, 8], [-45, 297], [237, 204], [234, 343], [374, -312], [30, -471], [106, -121], [301, 27], [93, -108], [137, -609], [317, -408], [181, -278], [291, -289], [369, -253], [-7, -362]], [[84713, 45326], [32, 139], [239, 133], [194, 20], [87, 74], [105, -74], [-102, -160], [-289, -258], [-233, -170]], [[32866, 56937], [160, 77], [58, -21], [-11, -440], [-232, -65], [-50, 53], [81, 163], [-6, 233]], [[52339, 72408], [302, 239], [195, -71], [-9, -299], [236, 217], [20, -113], [-139, -290], [-2, -273], [96, -147], [-36, -511], [-183, -297], [53, -322], [143, -10], [70, -281], [106, -92]], [[60041, 71744], [-102, 268], [105, 222], [-169, -51], [-233, 136], [-191, -340], [-421, -66], [-225, 317], [-300, 20], [-64, -245], [-192, -70], [-268, 314], [-303, -11], [-165, 588], [-203, 328], [135, 459], [-176, 283], [308, 565], [428, 23], [117, 449], [529, -78], [334, 383], [324, 167], [459, 13], [485, -417], [399, -228], [323, 91], [239, -53], [328, 309]], [[57776, 75399], [33, -228], [243, -190], [-51, -145], [-330, -33], [-118, -182], [-232, -319], [-87, 276], [3, 121]], [[83826, 64992], [-167, -947], [-119, -485], [-146, 499], [-32, 438], [163, 581], [223, 447], [127, -176], [-49, -357]], [[60889, 47817], [-128, -728], [16, -335], [178, -216], [8, -153], [-76, -357], [16, -180], [-18, -282], [97, -370], [115, -583], [101, -129]], [[59099, 45126], [-157, 177], [-177, 100], [-111, 99], [-116, 150]], [[58388, 46397], [-161, 331], [-55, 342]], [[58449, 49909], [98, 71], [304, -7], [566, 45]], [[60617, 78409], [-222, -48], [-185, -191], [-260, -31], [-239, -220], [16, -368], [136, -142], [284, 35], [-55, -210], [-304, -103], [-377, -342], [-154, 121], [61, 277], [-304, 173], [50, 113], [265, 197], [-80, 135], [-432, 149], [-19, 221], [-257, -73], [-103, -325], [-215, -437]], [[35174, 30629], [-121, -372], [-313, -328], [-205, 118], [-151, -63], [-256, 253], [-189, -19], [-169, 327]], [[6794, 61855], [-41, -99], [-69, 84], [8, 165], [-46, 216], [14, 65], [48, 97], [-19, 116], [16, 55], [21, -11], [107, -100], [49, -51], [45, -79], [71, -207], [-7, -33], [-108, -126], [-89, -92]], [[6645, 62777], [-94, -43], [-47, 125], [-32, 48], [-3, 37], [27, 50], [99, -56], [73, -90], [-23, -71]], [[6456, 63091], [-9, -63], [-149, 17], [21, 72], [137, -26]], [[6207, 63177], [-15, -34], [-19, 8], [-97, 21], [-35, 133], [-11, 24], [74, 82], [23, -38], [80, -196]], [[5737, 63567], [-33, -58], [-93, 107], [14, 43], [43, 58], [64, -12], [5, -138]], [[31350, 77248], [48, -194], [-296, -286], [-286, -204], [-293, -175], [-147, -351], [-47, -133], [-3, -313], [92, -313], [115, -15], [-29, 216], [83, -131], [-22, -169], [-188, -96], [-133, 11], [-205, -103], [-121, -29], [-162, -29], [-231, -171], [408, 111], [82, -112], [-389, -177], [-177, -1], [8, 72], [-84, -164], [82, -27], [-60, -424], [-203, -455], [-20, 152], [-61, 30], [-91, 148], [57, -318], [69, -105], [5, -223], [-89, -230], [-157, -472], [-25, 24], [86, 402], [-142, 225], [-33, 491], [-53, -255], [59, -375], [-183, 93], [191, -191], [12, -562], [79, -41], [29, -204], [39, -591], [-176, -439], [-288, -175], [-182, -346], [-139, -38], [-141, -217], [-39, -199], [-305, -383], [-157, -281], [-131, -351], [-43, -419], [50, -411], [92, -505], [124, -418], [1, -256], [132, -685], [-9, -398], [-12, -230], [-69, -361], [-83, -75], [-137, 72], [-44, 259], [-105, 136], [-148, 508], [-129, 452], [-42, 231], [57, 393], [-77, 325], [-217, 494], [-108, 90], [-281, -268], [-49, 30], [-135, 275], [-174, 147], [-314, -75], [-247, 66], [-212, -41], [-114, -92], [50, -157], [-5, -240], [59, -117], [-53, -77], [-103, 87], [-104, -112], [-202, 18], [-207, 312], [-242, -73], [-202, 137], [-173, -42], [-234, -138], [-253, -438], [-276, -255], [-152, -282], [-63, -266], [-3, -407], [14, -284], [52, -201]], [[17464, 69802], [-46, 302], [-180, 340], [-130, 71], [-30, 169], [-156, 30], [-100, 159], [-258, 59], [-71, 95], [-33, 324], [-270, 594], [-231, 821], [10, 137], [-123, 195], [-215, 495], [-38, 482], [-148, 323], [61, 489], [-10, 507], [-89, 453], [109, 557], [34, 536], [33, 536], [-50, 792], [-88, 506], [-80, 274], [33, 115], [402, -200], [148, -558], [69, 156], [-45, 484], [-94, 485]], [[7498, 84325], [-277, -225], [-142, 152], [-43, 277], [252, 210], [148, 90], [185, -40], [117, -183], [-240, -281]], [[4006, 85976], [-171, -92], [-182, 110], [-168, 161], [274, 101], [220, -54], [27, -226]], [[2297, 88264], [171, -113], [173, 61], [225, -156], [276, -79], [-23, -64], [-211, -125], [-211, 128], [-106, 107], [-245, -34], [-66, 52], [17, 223]], [[13740, 82958], [-153, 223], [-245, 188], [-78, 515], [-358, 478], [-150, 558], [-267, 38], [-441, 15], [-326, 170], [-574, 613], [-266, 112], [-486, 211], [-385, -51], [-546, 272], [-330, 252], [-309, -125], [58, -411], [-154, -38], [-321, -123], [-245, -199], [-308, -126], [-39, 348], [125, 580], [295, 182], [-76, 148], [-354, -329], [-190, -394], [-400, -420], [203, -287], [-262, -424], [-299, -248], [-278, -180], [-69, -261], [-434, -305], [-87, -278], [-325, -252], [-191, 45], [-259, -165], [-282, -201], [-231, -197], [-477, -169], [-43, 99], [304, 276], [271, 182], [296, 324], [345, 66], [137, 243], [385, 353], [62, 119], [205, 208], [48, 448], [141, 349], [-320, -179], [-90, 102], [-150, -215], [-181, 300], [-75, -212], [-104, 294], [-278, -236], [-170, 0], [-24, 352], [50, 216], [-179, 211], [-361, -113], [-235, 277], [-190, 142], [-1, 334], [-214, 252], [108, 340], [226, 330], [99, 303], [225, 43], [191, -94], [224, 285], [201, -51], [212, 183], [-52, 270], [-155, 106], [205, 228], [-170, -7], [-295, -128], [-85, -131], [-219, 131], [-392, -67], [-407, 142], [-117, 238], [-351, 343], [390, 247], [620, 289], [228, 0], [-38, -296], [586, 23], [-225, 366], [-342, 225], [-197, 296], [-267, 252], [-381, 187], [155, 309], [493, 19], [350, 270], [66, 287], [284, 281], [271, 68], [526, 262], [256, -40], [427, 315], [421, -124], [201, -266], [123, 114], [469, -35], [-16, -136], [425, -101], [283, 59], [585, -186], [534, -56], [214, -77], [370, 96], [421, -177], [302, -83]], [[30185, 57537], [-8, -139], [-163, -69], [91, -268], [-3, -309], [-123, -344], [105, -468], [120, 38], [62, 427], [-86, 208], [-14, 447], [346, 241], [-38, 278], [97, 186], [100, -415], [195, -9], [180, -330], [11, -195], [249, -6], [297, 61], [159, -264], [213, -74], [155, 185], [4, 149], [344, 35], [333, 9], [-236, -175], [95, -279], [222, -44], [210, -291], [45, -473], [144, 13], [109, -139]], [[80013, 63313], [-371, -505], [-231, -558], [-61, -410], [212, -623], [260, -772], [252, -365], [169, -475], [127, -1093], [-37, -1039], [-232, -389], [-318, -381], [-227, -492], [-346, -550], [-101, 378], [78, 401], [-206, 335]], [[96623, 40851], [-92, -78], [-93, 259], [10, 158], [175, -339]], [[96418, 41756], [45, -476], [-75, 74], [-58, -32], [-39, 163], [-6, 453], [133, -182]], [[64752, 60417], [-201, -158], [-54, -263], [-6, -201], [-277, -249], [-444, -276], [-249, -417], [-122, -33], [-83, 35], [-163, -245], [-177, -114], [-233, -30], [-70, -34], [-61, -156], [-73, -43], [-43, -150], [-137, 13], [-89, -80], [-192, 30], [-72, 345], [8, 323], [-46, 174], [-54, 437], [-80, 243], [56, 29], [-29, 270], [34, 114], [-12, 257]], [[58175, 37528], [113, -7], [134, -100], [94, 71], [148, -59]], [[59119, 34780], [-70, -430], [-32, -491], [-72, -267], [-190, -298], [-54, -86], [-118, -300], [-77, -303], [-158, -424], [-314, -609], [-196, -355], [-210, -269], [-290, -229], [-141, -31], [-36, -164], [-169, 88], [-138, -113], [-301, 114], [-168, -72], [-115, 31], [-286, -233], [-238, -94], [-171, -223], [-127, -14], [-117, 210], [-94, 11], [-120, 264], [-13, -82], [-37, 159], [2, 346], [-90, 396], [89, 108], [-7, 453], [-182, 553], [-139, 501], [-1, 1], [-199, 768]], [[58409, 41417], [-210, -81], [-159, -235], [-33, -205], [-100, -46], [-241, -486], [-154, -383], [-94, -13], [-90, 68], [-311, 65]]],
  "bbox": [-180, -85.60903777459767, 180, 83.64513000000001],
  "transform": {
    "scale": [0.0036000360003600037, 0.00169255860333201],
    "translate": [-180, -85.60903777459767]
  }
};
},{}],"index.js":[function(require,module,exports) {
"use strict";

var geo = _interopRequireWildcard(require("d3-geo"));

var select = _interopRequireWildcard(require("d3-selection"));

var _world110m = _interopRequireDefault(require("./assets/world-110m.json"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) { var desc = Object.defineProperty && Object.getOwnPropertyDescriptor ? Object.getOwnPropertyDescriptor(obj, key) : {}; if (desc.get || desc.set) { Object.defineProperty(newObj, key, desc); } else { newObj[key] = obj[key]; } } } } newObj.default = obj; return newObj; } }

console.log('hello world', geo, _world110m.default.objects.land); // const svg = document.getElementById('map')

var svg = select.select("#svg");
console.log(svg); // d3.json()
},{"d3-geo":"../node_modules/d3-geo/src/index.js","d3-selection":"../node_modules/d3-selection/src/index.js","./assets/world-110m.json":"assets/world-110m.json"}],"../node_modules/parcel-bundler/src/builtins/hmr-runtime.js":[function(require,module,exports) {
var global = arguments[3];
var OVERLAY_ID = '__parcel__error__overlay__';
var OldModule = module.bundle.Module;

function Module(moduleName) {
  OldModule.call(this, moduleName);
  this.hot = {
    data: module.bundle.hotData,
    _acceptCallbacks: [],
    _disposeCallbacks: [],
    accept: function (fn) {
      this._acceptCallbacks.push(fn || function () {});
    },
    dispose: function (fn) {
      this._disposeCallbacks.push(fn);
    }
  };
  module.bundle.hotData = null;
}

module.bundle.Module = Module;
var checkedAssets, assetsToAccept;
var parent = module.bundle.parent;

if ((!parent || !parent.isParcelRequire) && typeof WebSocket !== 'undefined') {
  var hostname = "" || location.hostname;
  var protocol = location.protocol === 'https:' ? 'wss' : 'ws';
  var ws = new WebSocket(protocol + '://' + hostname + ':' + "60095" + '/');

  ws.onmessage = function (event) {
    checkedAssets = {};
    assetsToAccept = [];
    var data = JSON.parse(event.data);

    if (data.type === 'update') {
      var handled = false;
      data.assets.forEach(function (asset) {
        if (!asset.isNew) {
          var didAccept = hmrAcceptCheck(global.parcelRequire, asset.id);

          if (didAccept) {
            handled = true;
          }
        }
      }); // Enable HMR for CSS by default.

      handled = handled || data.assets.every(function (asset) {
        return asset.type === 'css' && asset.generated.js;
      });

      if (handled) {
        console.clear();
        data.assets.forEach(function (asset) {
          hmrApply(global.parcelRequire, asset);
        });
        assetsToAccept.forEach(function (v) {
          hmrAcceptRun(v[0], v[1]);
        });
      } else {
        window.location.reload();
      }
    }

    if (data.type === 'reload') {
      ws.close();

      ws.onclose = function () {
        location.reload();
      };
    }

    if (data.type === 'error-resolved') {
      console.log('[parcel] ✨ Error resolved');
      removeErrorOverlay();
    }

    if (data.type === 'error') {
      console.error('[parcel] 🚨  ' + data.error.message + '\n' + data.error.stack);
      removeErrorOverlay();
      var overlay = createErrorOverlay(data);
      document.body.appendChild(overlay);
    }
  };
}

function removeErrorOverlay() {
  var overlay = document.getElementById(OVERLAY_ID);

  if (overlay) {
    overlay.remove();
  }
}

function createErrorOverlay(data) {
  var overlay = document.createElement('div');
  overlay.id = OVERLAY_ID; // html encode message and stack trace

  var message = document.createElement('div');
  var stackTrace = document.createElement('pre');
  message.innerText = data.error.message;
  stackTrace.innerText = data.error.stack;
  overlay.innerHTML = '<div style="background: black; font-size: 16px; color: white; position: fixed; height: 100%; width: 100%; top: 0px; left: 0px; padding: 30px; opacity: 0.85; font-family: Menlo, Consolas, monospace; z-index: 9999;">' + '<span style="background: red; padding: 2px 4px; border-radius: 2px;">ERROR</span>' + '<span style="top: 2px; margin-left: 5px; position: relative;">🚨</span>' + '<div style="font-size: 18px; font-weight: bold; margin-top: 20px;">' + message.innerHTML + '</div>' + '<pre>' + stackTrace.innerHTML + '</pre>' + '</div>';
  return overlay;
}

function getParents(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return [];
  }

  var parents = [];
  var k, d, dep;

  for (k in modules) {
    for (d in modules[k][1]) {
      dep = modules[k][1][d];

      if (dep === id || Array.isArray(dep) && dep[dep.length - 1] === id) {
        parents.push(k);
      }
    }
  }

  if (bundle.parent) {
    parents = parents.concat(getParents(bundle.parent, id));
  }

  return parents;
}

function hmrApply(bundle, asset) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (modules[asset.id] || !bundle.parent) {
    var fn = new Function('require', 'module', 'exports', asset.generated.js);
    asset.isNew = !modules[asset.id];
    modules[asset.id] = [fn, asset.deps];
  } else if (bundle.parent) {
    hmrApply(bundle.parent, asset);
  }
}

function hmrAcceptCheck(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (!modules[id] && bundle.parent) {
    return hmrAcceptCheck(bundle.parent, id);
  }

  if (checkedAssets[id]) {
    return;
  }

  checkedAssets[id] = true;
  var cached = bundle.cache[id];
  assetsToAccept.push([bundle, id]);

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    return true;
  }

  return getParents(global.parcelRequire, id).some(function (id) {
    return hmrAcceptCheck(global.parcelRequire, id);
  });
}

function hmrAcceptRun(bundle, id) {
  var cached = bundle.cache[id];
  bundle.hotData = {};

  if (cached) {
    cached.hot.data = bundle.hotData;
  }

  if (cached && cached.hot && cached.hot._disposeCallbacks.length) {
    cached.hot._disposeCallbacks.forEach(function (cb) {
      cb(bundle.hotData);
    });
  }

  delete bundle.cache[id];
  bundle(id);
  cached = bundle.cache[id];

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    cached.hot._acceptCallbacks.forEach(function (cb) {
      cb();
    });

    return true;
  }
}
},{}]},{},["../node_modules/parcel-bundler/src/builtins/hmr-runtime.js","index.js"], null)
//# sourceMappingURL=/src.e31bb0bc.js.map