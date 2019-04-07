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
},{"./create":"../node_modules/d3-selection/src/create.js","./creator":"../node_modules/d3-selection/src/creator.js","./local":"../node_modules/d3-selection/src/local.js","./matcher":"../node_modules/d3-selection/src/matcher.js","./mouse":"../node_modules/d3-selection/src/mouse.js","./namespace":"../node_modules/d3-selection/src/namespace.js","./namespaces":"../node_modules/d3-selection/src/namespaces.js","./point":"../node_modules/d3-selection/src/point.js","./select":"../node_modules/d3-selection/src/select.js","./selectAll":"../node_modules/d3-selection/src/selectAll.js","./selection/index":"../node_modules/d3-selection/src/selection/index.js","./selector":"../node_modules/d3-selection/src/selector.js","./selectorAll":"../node_modules/d3-selection/src/selectorAll.js","./selection/style":"../node_modules/d3-selection/src/selection/style.js","./touch":"../node_modules/d3-selection/src/touch.js","./touches":"../node_modules/d3-selection/src/touches.js","./window":"../node_modules/d3-selection/src/window.js","./selection/on":"../node_modules/d3-selection/src/selection/on.js"}],"assets/ne_110m_land.json":[function(require,module,exports) {
module.exports = {
  "type": "FeatureCollection",
  "features": [{
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-59.572094692611529, -80.040178725096297], [-59.865849371974633, -80.549656671061868], [-60.159655727770115, -81.000326837079314], [-62.255393439367083, -80.863177585776668], [-64.488125372969847, -80.921933689292572], [-65.74166642928995, -80.588827406739142], [-65.74166642928995, -80.549656671061868], [-66.290030890555045, -80.25577280061799], [-64.037687750897646, -80.294943536295278], [-61.883245612217053, -80.392870375488286], [-61.138975796133366, -79.9813709451481], [-60.61011918805832, -79.628679294756125], [-59.572094692611529, -80.040178725096297]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-159.208183560197654, -79.497059421708727], [-161.127601284814773, -79.634208673011301], [-162.439846768218388, -79.281465346187034], [-163.027407803376974, -78.928773695794973], [-163.066604377270437, -78.869965915846848], [-163.712895677728739, -78.595667413241529], [-163.712895677728739, -78.595666605797291], [-163.105800951163815, -78.223337911134394], [-161.245113491846411, -78.380175883140168], [-160.246208055644502, -78.693645121422676], [-159.482404548154477, -79.046337579258989], [-159.208183560197654, -79.497059421708727]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-45.154757656421026, -78.047069600586738], [-43.920827806155671, -78.478102722333261], [-43.489949713706096, -79.085559991368541], [-43.372437506674459, -79.516644789547371], [-43.333266770997113, -80.026122735512928], [-44.880536668464288, -80.339643650227714], [-46.506173875501929, -80.594356784994318], [-48.386420864441874, -80.829484551922363], [-50.482106899606464, -81.025441583173134], [-52.851988084511703, -80.966685479657315], [-54.164259406131521, -80.633527520671592], [-53.987991095583965, -80.222028090331406], [-51.8531343247422, -79.947729587726087], [-50.991326463410587, -79.614623305172756], [-50.364594692574656, -79.183486830561648], [-49.914131232286451, -78.811209004886706], [-49.306958991073117, -78.458569030926924], [-48.660616014182438, -78.04701792415446], [-48.660616014182438, -78.047018731598698], [-48.151396450378314, -78.047069600586738], [-46.66285681821094, -77.831475525065045], [-45.154757656421026, -78.047069600586738]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-121.211511393857137, -73.500990499006051], [-119.918851278292053, -73.65772511814734], [-118.724143032692012, -73.481353454735213], [-119.292118700011955, -73.834096781559481], [-120.232217163709976, -74.08880991632617], [-121.622829956684257, -74.010468444971707], [-122.621734585441899, -73.657777602023884], [-122.621735392886237, -73.657776794579632], [-122.406244670229114, -73.324618835593924], [-121.211511393857137, -73.500990499006051]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-125.559566406895328, -73.481353454735213], [-124.031881877266869, -73.873267517236755], [-124.619468750641545, -73.834096781559481], [-125.912180542638907, -73.736118265934095], [-127.283129645681925, -73.461768894340821], [-127.283130453126248, -73.461768086896569], [-126.558471843097308, -73.246225687807168], [-125.559566406895328, -73.481353454735213]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-98.981549648823915, -71.933334248999785], [-97.884743211645116, -72.070535176734751], [-96.787936774466232, -71.952971293270721], [-96.200349901091471, -72.521205342752197], [-96.983764614636243, -72.442863871397648], [-98.19808325884685, -72.482034607074922], [-99.43201310911212, -72.442863871397648], [-100.783455166409212, -72.501619974913552], [-101.801868455801355, -72.305662943662782], [-102.3307250638764, -71.894164320766862], [-102.3307250638764, -71.89416351332261], [-101.703967454824451, -71.717791849910384], [-100.4309185453141, -71.854992777645336], [-98.981549648823915, -71.933334248999785]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-68.451345994730332, -70.955822855766741], [-68.333833787698779, -71.4064930217842], [-68.510127936462339, -71.798407084285742], [-68.784297247986984, -72.170735778948639], [-69.959470994736421, -72.307885030251299], [-71.075888637970053, -72.50384206150207], [-72.388134121373781, -72.484256693663525], [-71.898499925408203, -72.092342631161898], [-73.073621995725546, -72.229491882464544], [-74.190039638959064, -72.366692810199496], [-74.953894822881495, -72.072757263323254], [-75.01262508818121, -71.661257832983068], [-73.915818651002326, -71.269344577925779], [-73.915818651002326, -71.269343770481527], [-73.230330776650646, -71.151779887017511], [-72.074716559523466, -71.190950622694785], [-71.780961880160362, -70.681472676729214], [-71.722179938428354, -70.30919565849851], [-71.741791144483187, -69.505782165656797], [-71.173815477163146, -69.035474955368414], [-70.253251512315813, -68.87874033622721], [-69.72444658067306, -69.251017354457815], [-69.489422166609586, -69.623346049120812], [-69.058518235943836, -70.074016215138187], [-68.725541144471066, -70.505152689749281], [-68.451345994730332, -70.955822855766741]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-58.614142829001025, -64.152467130133232], [-59.045072597882893, -64.368009529222633], [-59.78934241396658, -64.21122323364915], [-60.611927863188612, -64.309201749274536], [-61.297415737540376, -64.544329516202481], [-62.022100185785433, -64.799094327401463], [-62.511760219967101, -65.093029874277619], [-62.648857794837369, -65.484942321890657], [-62.590127529537654, -65.857219340121361], [-62.120078701410904, -66.190325622674806], [-62.805566575762583, -66.425505066035043], [-63.745690070232456, -66.503846537389592], [-64.294106207929957, -66.8370044963753], [-64.881693081304718, -67.150473734657723], [-65.508424852140479, -67.581610209268831], [-65.665081956633401, -67.953887227499536], [-65.31254533553809, -68.365334981407415], [-64.783714565679333, -68.678907572554508], [-63.961103278241097, -68.91398366305026], [-63.197299770750988, -69.227556254197339], [-62.78595536970775, -69.619418640266588], [-62.570516323482934, -69.991747334929585], [-62.276735805903627, -70.383661397431126], [-61.806661139560589, -70.716767679984571], [-61.512906460197399, -71.089044698215176], [-61.375808885327132, -72.010073750953211], [-61.081976691315447, -72.382350769183915], [-61.003661058177187, -72.774264831685457], [-60.69026933454316, -73.166178894186999], [-60.827366909413456, -73.6952422079912], [-61.375808885327132, -74.106741638331471], [-61.963369920485604, -74.439847920884802], [-63.295200771727963, -74.576997172187475], [-63.745690070232456, -74.929740499011814], [-64.352836473229672, -75.262846781565258], [-65.860987311451879, -75.635123799795878], [-67.192818162694039, -75.791910095369445], [-68.446281704365845, -76.00745249445886], [-69.797723761662837, -76.222994893548275], [-70.600723843046353, -76.634494323888447], [-72.206775682245365, -76.673665059565636], [-73.969536302369761, -76.634494323888447], [-75.555976935513939, -76.712887471675202], [-77.240370246067698, -76.712887471675202], [-76.926978522433672, -77.104801534176744], [-75.399293992805127, -77.2810698447243], [-74.282876349571495, -77.555420023761911], [-73.656118740519361, -77.908111674153872], [-74.772536383753078, -78.221632588868772], [-76.496100429984011, -78.123654073243287], [-77.925858120419377, -78.378418884442254], [-77.984665900367474, -78.78991831478244], [-78.02378495961247, -79.181833184728305], [-76.848637051079123, -79.514939467281664], [-76.633223843070482, -79.887216485512369], [-75.360097418911749, -80.259545180175266], [-73.244851854124619, -80.416331475748848], [-71.442946336539222, -80.690629978354082], [-70.013162807887767, -81.004150893068868], [-68.191646084247651, -81.317671807783668], [-65.704278530526665, -81.47445810335725], [-63.256030036050873, -81.748756605962484], [-61.552025519442481, -82.04269215283864], [-59.691415574773458, -82.375850111824349], [-58.712121344626411, -82.846105645680439], [-58.222487148660832, -83.21843434034335], [-57.008116828018018, -82.865691013519083], [-55.362894253141633, -82.571755466642841], [-53.619770677288244, -82.258234551928041], [-51.543644171746109, -82.003521417161451], [-49.761349860215546, -81.72917123812384], [-47.273930630062381, -81.709585870285295], [-44.825707973802594, -81.846735121587869], [-42.808363409992467, -82.081914564948192], [-42.162020433101873, -81.650829766769391], [-40.771433478343596, -81.356894219893235], [-38.244817674296968, -81.337308852054591], [-36.266669684380304, -81.121714776532897], [-34.38639685722427, -80.906172377443482], [-32.31029618989831, -80.769023126140922], [-30.097097947701997, -80.592651462728696], [-28.549802212018704, -80.337938327962092], [-29.254901292425131, -79.985195001137754], [-29.685805223090966, -79.632503350745679], [-29.685805223090966, -79.26022633251506], [-31.624808315546545, -79.299397068192249], [-33.681323615034046, -79.456131687333539], [-35.639912075328283, -79.456131687333539], [-35.914107225069017, -79.083854669102919], [-35.777009650198636, -78.339248148764995], [-35.326546189910431, -78.123654073243287], [-33.896762661258805, -77.888526306315327], [-32.212369350705217, -77.653450215819575], [-30.99805070649461, -77.359514668943419], [-29.783732062284059, -77.065579122067177], [-28.88277930349139, -76.673665059565636], [-27.511751878355653, -76.497345072585787], [-26.16033565927475, -76.360144144850835], [-25.474821946706868, -76.281802673496372], [-23.927552049239779, -76.24258026138682], [-22.45859778491095, -76.105431010084246], [-21.224693772861798, -75.90947397883339], [-20.010375128651162, -75.67434621190543], [-18.913542853256189, -75.439218444977385], [-17.522981736714115, -75.125697530262585], [-16.641588507543958, -74.792539571276876], [-15.701490851290174, -74.49860402440072], [-15.407710333710867, -74.106741638331471], [-16.465320196996402, -73.871613871403412], [-16.112783575901176, -73.46011444106324], [-15.446855231171952, -73.146541849916161], [-14.408804897508986, -72.950584818665305], [-13.311972622113984, -72.715457051737346], [-12.293507656289563, -72.401936137022545], [-11.510067104528616, -72.010073750953211], [-11.020432908563123, -71.539766540664914], [-10.295774298534184, -71.265416361627302], [-9.101015183946089, -71.324224141575513], [-8.611380987980596, -71.657330424128943], [-7.416621873392444, -71.696501159806132], [-7.377451137715155, -71.324224141575513], [-6.868231573911117, -70.932310079073972], [-5.790984666354689, -71.030288594699357], [-5.53637488445267, -71.402617289362269], [-4.341667446296782, -71.461373392878173], [-3.048981492515622, -71.285053405898239], [-1.795492112627812, -71.167437846001917], [-0.65948910155555, -71.226245625950128], [-0.228636847322065, -71.637745056290299], [0.868195428072994, -71.304638773736883], [1.886686232113533, -71.128267110324742], [3.022637566753417, -70.991117859022083], [4.139055209987049, -70.853916931287131], [5.157546014027673, -70.618789164359185], [6.273911980828984, -70.462054545217882], [7.135719842160512, -70.246512146128481], [7.742866245157728, -69.893768819304114], [8.487110223025326, -70.148533630503096], [9.525134718472202, -70.011332702768129], [10.249845004933434, -70.481639913056512], [10.81782067225339, -70.834331563448501], [11.953823683325652, -70.63837453219773], [12.404287143613942, -70.246512146128481], [13.422777947654453, -69.972161967090955], [14.734997592842006, -70.030918070606774], [15.126756626046671, -70.40324676526977], [15.949342075268703, -70.030918070606774], [17.026588982825075, -69.913354187142744], [18.201711053142304, -69.874183451465569], [19.259372592860103, -69.893768819304114], [20.375738559661443, -70.011332702768129], [21.452985467217815, -70.07014048271634], [21.92303429534465, -70.40324676526977], [22.569403110451418, -70.697182312145927], [23.666183709414128, -70.5208106487338], [24.841357456163649, -70.481639913056512], [25.977308790803647, -70.481639913056512], [27.093726434037251, -70.462054545217882], [28.092580193806867, -70.32485361748293], [29.150241733524666, -70.2072897340189], [30.031583286262531, -69.932939554981388], [30.971732618948607, -69.756619568001554], [31.990171746556854, -69.658641052376169], [32.754052768695345, -69.384290873338557], [33.302443068176643, -68.835642191695712], [33.870418735496713, -68.502587585574659], [34.90849490737574, -68.659270528283571], [35.300202264148169, -69.012013855107938], [36.162010125479782, -69.247141622035983], [37.20003462092663, -69.168748474249142], [37.905107863116967, -69.521440124641202], [38.6494035174168, -69.776204935840099], [39.667894321457339, -69.541077168912139], [40.020430942552451, -69.109940694301031], [40.921357863129089, -68.933620707321182], [41.959434035008115, -68.600514424767766], [42.938702426939216, -68.4633134970328], [44.113876173688624, -68.267408142214322], [44.897290887233481, -68.051865743124921], [45.719928012887834, -67.816737976196876], [46.503342726432635, -67.601195577107461], [47.443440382686362, -67.718759460571491], [48.344418979695178, -67.36606781017943], [48.990736118369597, -67.091717631141904], [49.930885451055673, -67.111302998980534], [50.753470900277676, -66.876175232052503], [50.949324578663919, -66.523483581660429], [51.791547072156874, -66.249133402622903], [52.614132521378934, -66.053176371372132], [53.613037957580872, -65.89639007579855], [54.53355024599594, -65.818048604444101], [55.414943475166211, -65.876804707959906], [56.355041131419995, -65.974783223585391], [57.158092889235689, -66.249133402622903], [57.255968051996462, -66.680218200801733], [58.137361281166619, -67.013324483355149], [58.744507684163835, -67.287674662392675], [59.939318475184308, -67.40523854585669], [60.605220981697329, -67.679588724894217], [61.427806430919389, -67.953887227499536], [62.387489455011831, -68.012695007447647], [63.190489536395233, -67.816737976196876], [64.052349074159139, -67.40523854585669], [64.992446730412809, -67.620729268513813], [65.971715122343852, -67.738344828410135], [66.911864455029871, -67.85590871187415], [67.891132846960801, -67.934301859660906], [68.890038283162795, -67.934301859660906], [69.712623732384856, -68.972791442998371], [69.673452996707539, -69.227556254197339], [69.555940789675816, -69.678226420214713], [68.596257765583488, -69.932939554981388], [67.812739699174045, -70.305268249644286], [67.949888950476719, -70.697182312145927], [69.066306593710323, -70.67754526787499], [68.929157342407706, -71.069459330376532], [68.419989455035847, -71.441788025039529], [67.949888950476719, -71.853287455379714], [68.713769972615211, -72.166808370094515], [69.869306675093782, -72.2647868857199], [71.024895054004645, -72.088415222307759], [71.573285353486057, -71.696501159806132], [71.906288283174803, -71.324224141575513], [72.454626906223922, -71.010703226860727], [73.081410353492089, -70.716767679984571], [73.336020135394193, -70.364024353160204], [73.864876743469125, -69.874183451465569], [74.491556837872707, -69.776204935840099], [75.627559848944969, -69.73703420016291], [76.62646528514685, -69.619418640266588], [77.644904412755096, -69.462684021125398], [78.13453860872059, -69.070769958623757], [78.428370802732132, -68.698441263960859], [79.113858677083925, -68.326215922162447], [80.093127069014855, -68.071502787395843], [80.935349562507866, -67.875545756145073], [81.483791538421514, -67.542387797159364], [82.05176720574147, -67.36606781017943], [82.776425815770409, -67.209281514605919], [83.77533125197229, -67.307260030231305], [84.676206496116635, -67.209281514605919], [85.655526564479857, -67.091717631141904], [86.752358839874859, -67.150473734657723], [87.477017449903798, -66.876175232052503], [87.986288690140242, -66.20991099051335], [88.358410679073955, -66.484261169550962], [88.828407830768498, -66.954568379839245], [89.670630324261509, -67.150473734657723], [90.630365024786187, -67.228866882444564], [91.590099725310807, -67.111302998980534], [92.608538852919168, -67.18969614676729], [93.548636509172951, -67.209281514605919], [94.175419956441004, -67.111302998980534], [95.017590773501666, -67.17011077892866], [95.781471795640158, -67.38565317801806], [96.682398716216682, -67.248503926715486], [97.759645623773025, -67.248503926715486], [98.680209588620443, -67.111302998980534], [99.718182407635055, -67.248503926715486], [100.384188267012661, -66.915345967729763], [100.89335615438452, -66.582239685176347], [101.578895705168549, -66.307889506138736], [102.832410923272533, -65.563283793245205], [103.47867638551466, -65.700484720980072], [104.242557407653152, -65.974783223585391], [104.90845991416623, -66.327526550409658], [106.181560500108816, -66.934931335568308], [107.160880568472209, -66.954568379839245], [108.08139285688722, -66.954568379839245], [109.158639764443677, -66.8370044963753], [110.235834995567728, -66.699803568640363], [111.058472121222138, -66.425505066035043], [111.743959995573817, -66.131569519158887], [112.860377638807364, -66.092347107049321], [113.604673293107368, -65.876804707959906], [114.388088006651998, -66.072761739210776], [114.897307570456206, -66.386282653925576], [115.60238081264643, -66.699803568640363], [116.699161411609253, -66.660632832963088], [117.384700962393225, -66.915345967729763], [118.579460076981292, -67.17011077892866], [119.832923618652984, -67.268089294554045], [120.870999790532181, -67.18969614676729], [121.654414504076925, -66.876175232052503], [122.320368687022352, -66.562654317337703], [123.221295607598876, -66.484261169550962], [124.122274204607635, -66.621462097285814], [125.160247023622304, -66.719388936478992], [126.100396356308323, -66.562654317337703], [127.001426629749375, -66.562654317337703], [127.88276818248724, -66.660632832963088], [128.803280470902422, -66.758611348588474], [129.704259067911181, -66.582239685176347], [130.781454299035346, -66.425505066035043], [131.799945103075942, -66.386282653925576], [132.935896437715854, -66.386282653925576], [133.856460402563272, -66.288304138300191], [134.757387323139795, -66.209962666945643], [135.031582472880444, -65.720070088818716], [135.070753208557704, -65.30857065847853], [135.697484979393579, -65.58286916108375], [135.873804966373399, -66.033591003533502], [136.206704543197674, -66.445090433873673], [136.618048944240996, -66.778196716427118], [137.460271437733951, -66.954568379839245], [138.59622277237392, -66.895760599891133], [139.908442417561361, -66.876175232052503], [140.809421014570177, -66.817367452104378], [142.12169233619008, -66.817367452104378], [143.061841668876156, -66.797782084265748], [144.374061314063709, -66.8370044963753], [145.490427280865021, -66.915345967729763], [146.195552199487594, -67.228866882444564], [145.999698521101408, -67.601195577107461], [146.646067336208233, -67.895131123983617], [147.723262567332284, -68.130258890911676], [148.839628534133595, -68.38502370211063], [150.132314487914783, -68.561292012658186], [151.483704868779654, -68.718129984664074], [152.502247349252428, -68.874812927372986], [153.638198683892455, -68.894501648076201], [154.284567498999166, -68.561292012658186], [155.165857375304739, -68.835642191695712], [155.929790073875523, -69.149214782842876], [156.811131626613388, -69.384290873338557], [158.025527785472406, -69.482269388963942], [159.181012811518741, -69.599833272427958], [159.670698683916584, -69.991747334929585], [160.806650018556496, -70.226875101857544], [161.570479364262695, -70.579618428681897], [162.686897007496242, -70.736353047823116], [163.84243370997487, -70.716767679984571], [164.919680617531213, -70.775523783500375], [166.114439732119394, -70.755938415661745], [167.309095493842875, -70.834331563448501], [168.425616489941063, -70.971480814751146], [169.463589308955676, -71.206660258111498], [170.50166548083493, -71.402617289362269], [171.206790399457446, -71.696501159806132], [171.089226515993488, -72.088415222307759], [170.560421584350735, -72.441158549132027], [170.109958124062445, -72.891828715149472], [169.75736982653504, -73.244520365541547], [169.287320998408092, -73.656019795881718], [167.975101353220595, -73.812806091455229], [167.387488641629631, -74.16549774184729], [166.094802687848443, -74.381040140936705], [165.644390903992445, -74.772954203438246], [164.958851353208473, -75.145282898101229], [164.234192743179534, -75.458803812816029], [163.822796665703919, -75.870303243156201], [163.568238560234221, -76.24258026138682], [163.470260044608807, -76.693302103836558], [163.489897088879815, -77.065579122067177], [164.057872756199771, -77.457441508136526], [164.273363478856908, -77.829770202799423], [164.743463983416035, -78.182513529623776], [166.604125604517122, -78.319611104494157], [166.995781284857316, -78.750747579105166], [165.193875767272033, -78.907483005690708], [163.666217075859578, -79.123025404780122], [161.766384719081174, -79.162247816889675], [160.924162225588219, -79.730481866371065], [160.747893915040578, -80.200737400227155], [160.31696414615871, -80.573066094890066], [159.78821089094825, -80.945394789553049], [161.120015903974405, -81.278501072106479], [161.629287144210906, -81.690000502446651], [162.490991652677764, -82.06227752067727], [163.705336135104545, -82.395435479662993], [165.095948928078855, -82.708956394377793], [166.604125604517122, -83.022477309092579], [168.895665318068012, -83.335998223807366], [169.404781529007636, -83.825890801934392], [172.283933954149319, -84.041433201023793], [172.477048781624006, -84.117914320815686], [173.224083286835395, -84.413710219254412], [175.985671828513063, -84.158997084487737], [178.277211542064066, -84.472517999202523], [180.000000000000142, -84.71338], [180.000000000000142, -90.0], [-180.0, -90.0], [-180.0, -84.71338], [-179.942499356179042, -84.721443373552489], [-179.058677334691197, -84.139411716649178], [-177.256771817105829, -84.452932631363893], [-177.140806673265871, -84.417941227148319], [-176.861992942389094, -84.333811995377147], [-176.523951560551779, -84.2318107924755], [-176.230303463830239, -84.143203474868585], [-176.084672818077706, -84.099259128758334], [-175.934100613487232, -84.101591027765551], [-175.829882168662635, -84.117914320815686], [-174.382502814815695, -84.534323012223652], [-173.116559414745467, -84.117914320815686], [-172.889105598012804, -84.061018568862352], [-169.951222907571349, -83.884646905450211], [-168.999988980158747, -84.117914320815686], [-168.530198534193346, -84.237390232274564], [-167.022099372403432, -84.570496514827909], [-164.182143521155155, -84.825209649594598], [-161.929774543281468, -85.138730564309384], [-158.071379564424859, -85.373910007669721], [-155.192252977499209, -85.099559828632195], [-150.942098965438021, -85.295516859882966], [-148.533072883071611, -85.609037774597766], [-145.888918226332976, -85.31510222772161], [-143.107718478600447, -85.040752048683999], [-142.892279432375631, -84.570496514827909], [-146.829068366463304, -84.531274102718442], [-150.060731574483952, -84.296146335790382], [-150.902928229760846, -83.904232273288841], [-153.586201138300197, -83.68868987419944], [-153.409906989536466, -83.23801970818198], [-153.037759162386521, -82.826520277841809], [-152.665637173452751, -82.454191583178897], [-152.861516690055055, -82.04269215283864], [-154.526298794553981, -81.768393650233406], [-155.290179816692387, -81.415650323409054], [-156.83744971415959, -81.102129408694253], [-154.408786587522229, -81.160937188642464], [-152.097661506132823, -81.004150893068868], [-150.648292609642624, -81.337308852054591], [-148.865998298112061, -81.043373305178434], [-147.220749885019501, -80.671044610515452], [-146.417748996191847, -80.337938327962092], [-146.770286424731296, -79.926438897621836], [-148.062946540296281, -79.652088718584309], [-149.53190080462511, -79.358204848140446], [-151.588416104112525, -79.299397068192249], [-153.390321621697808, -79.162247816889675], [-155.329376390585765, -79.064269301264289], [-155.97566769104418, -78.691939799157055], [-157.268301968393075, -78.378418884442254], [-158.051768358370111, -78.025675557617902], [-158.365134243788049, -76.889207458655051], [-157.875474209606466, -76.987237650712714], [-156.97457312724606, -77.300758565427515], [-155.329376390585765, -77.202728373369837], [-153.742832404576831, -77.065579122067177], [-152.920246955354685, -77.496663920245993], [-151.333780483994303, -77.3987370810529], [-150.00194963275186, -77.183143005531207], [-148.748486091080338, -76.908844502925973], [-147.612483080008076, -76.575738220372529], [-146.104408948990113, -76.477759704747143], [-146.14352800823491, -76.105431010084246], [-146.496091274990533, -75.733153991853541], [-146.202309949967002, -75.380410665029189], [-144.909623996185843, -75.204039001617048], [-144.322037122811082, -75.53719696060277], [-142.794352593182509, -75.341239929352], [-141.638764214271703, -75.086475118153032], [-140.209006523836251, -75.066889750314488], [-138.857590304755377, -74.968911234689017], [-137.506199923890563, -74.733783467761043], [-136.428901339901927, -74.518241068671642], [-135.214582695691291, -74.302698669582242], [-134.431193820362523, -74.36145477309806], [-133.745654269578665, -74.439847920884802], [-132.257167928732059, -74.302698669582242], [-130.925311239273611, -74.47901865656209], [-129.554283814137875, -74.459433288723446], [-128.242038330734147, -74.322284037420786], [-126.890622111653258, -74.420262553046257], [-125.402082479485884, -74.518241068671642], [-124.011495524727707, -74.47901865656209], [-122.562152466453711, -74.49860402440072], [-121.073612834286237, -74.518241068671642], [-119.702559570934412, -74.47901865656209], [-118.684145474098031, -74.185083109685934], [-117.469800991671306, -74.02834849054463], [-116.216311611783496, -74.243890889634031], [-115.021552497195415, -74.067519226221904], [-113.944331427855161, -73.714827575829844], [-113.297988450964482, -74.02834849054463], [-112.94545182986937, -74.381040140936705], [-112.299083014762687, -74.714198099922413], [-111.26105851931581, -74.420262553046257], [-110.066325242943734, -74.792539571276876], [-108.71490902386283, -74.910103454740891], [-107.559346483168127, -75.184453633778418], [-106.149148322355117, -75.125697530262585], [-104.876073574628762, -74.949325866850458], [-103.367948574622758, -74.988496602527647], [-102.016506517325666, -75.125697530262585], [-100.645530768622308, -75.302017517242433], [-100.116699998763366, -74.870932719063632], [-100.76304297565396, -74.537826436510272], [-101.252703009835628, -74.185083109685934], [-102.545337287184608, -74.106741638331471], [-103.113312954504551, -73.734412943668474], [-103.328752000729381, -73.362084249005562], [-103.681288621824493, -72.617530212544239], [-102.917485114334454, -72.754679463846827], [-101.605239630930825, -72.813435567362731], [-100.312527838933462, -72.754679463846827], [-99.137379930400115, -72.911414082988117], [-98.118889126359591, -73.205349629864273], [-97.68803687212602, -73.558041280256333], [-96.336594814828942, -73.616849060204459], [-95.043960537480046, -73.47969980890187], [-93.672907274128221, -73.283742777651014], [-92.439003262079041, -73.166178894186999], [-91.420564134470709, -73.401306661115044], [-90.088733283228351, -73.322913513328302], [-89.226951260113026, -72.558722432596056], [-88.423951178729624, -73.009392598613502], [-87.268336961602614, -73.185764262025629], [-86.014821743498629, -73.087785746400243], [-85.192236294276483, -73.47969980890187], [-83.879990810872954, -73.518870544579073], [-82.665646328446229, -73.636434428043088], [-81.470913052074138, -73.851976827132489], [-80.687446662096988, -73.47969980890187], [-80.295790981756994, -73.126956482077532], [-79.296885545555, -73.518870544579073], [-77.925858120419377, -73.420892028953688], [-76.907367316378753, -73.636434428043088], [-76.221879442027074, -73.969540710596519], [-74.890048590784801, -73.871613871403412], [-73.852024095337924, -73.656019795881718], [-72.833533291297499, -73.401306661115044], [-71.619214647086778, -73.264157409812469], [-70.209042324490071, -73.146541849916161], [-68.935915900331338, -73.009392598613502], [-67.956621670184262, -72.793850199524087], [-67.369060635025505, -72.480329284809301], [-67.134036220962031, -72.049244486630485], [-67.251548427993868, -71.637745056290299], [-67.56494015162798, -71.245830993788758], [-67.917476772723091, -70.853916931287131], [-68.230842658140915, -70.462054545217882], [-68.485452440043019, -70.109311218393515], [-68.544208543558852, -69.717397155891973], [-68.446281704365845, -69.325534769822724], [-67.97623287623901, -68.953206075159727], [-67.584499681250435, -68.541706644819556], [-67.427842576757513, -68.149844258750306], [-67.623670416927609, -67.718759460571491], [-67.741182623959247, -67.326845398069949], [-67.251548427993868, -66.876175232052503], [-66.703183966728659, -66.582239685176347], [-66.056815151621976, -66.209962666945643], [-65.371327277270211, -65.89639007579855], [-64.568275519454403, -65.602506205354686], [-64.176542324465828, -65.17142302206436], [-63.62815202498453, -64.897072843026848], [-63.001394415932481, -64.642308031827866], [-62.041685553624063, -64.583551928312048], [-61.414927944572014, -64.270031013597247], [-60.709854702381705, -64.074073982346476], [-59.887269253159758, -63.956510098882461], [-59.162584804914701, -63.701745287683487], [-58.594557461162282, -63.388224372968693], [-57.811142747617509, -63.270660489504671], [-57.223581712459037, -63.525425300703645], [-57.595729539608783, -63.858531583257076], [-58.614142829001025, -64.152467130133232]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-67.75, -53.85], [-66.45, -54.45], [-65.05, -54.699999999999882], [-65.5, -55.2], [-66.45, -55.25], [-66.95992, -54.896809999999888], [-67.291029999999893, -55.30124], [-68.14863, -55.611829999999877], [-69.2321, -55.49906], [-69.95809, -55.19843], [-71.00568, -55.05383], [-72.2639, -54.49514], [-73.2852, -53.957519999999896], [-74.66253, -52.83749], [-73.8381, -53.047429999999899], [-72.43418, -53.7154], [-71.10773, -54.07433], [-70.59178, -53.61583], [-70.26748, -52.93123], [-69.34565, -52.5183], [-68.63411, -52.63625], [-68.63401022758319, -52.636370458874374], [-68.25, -53.099999999999881], [-67.75, -53.85]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-58.55, -51.1], [-57.75, -51.55], [-58.05, -51.9], [-59.4, -52.2], [-59.85, -51.85], [-60.7, -52.3], [-61.2, -51.85], [-60.0, -51.25], [-59.15, -51.5], [-58.55, -51.1]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[70.28, -49.71], [68.745, -49.775], [68.72, -49.2425], [68.8675, -48.83], [68.935, -48.625], [69.580000000000126, -48.94], [70.525, -49.065], [70.56, -49.255], [70.28, -49.71]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[145.397978143494896, -40.79254851660599], [146.364120721623664, -41.137695407883342], [146.908583612250851, -41.000546156580782], [147.689259474884153, -40.808258152022688], [148.289067824496072, -40.87543751400213], [148.359864536735898, -42.062445163746446], [148.017301467073111, -42.407023614268624], [147.91405195535387, -43.211522312188507], [147.564564243764011, -42.937688897473876], [146.870343052354968, -43.634597263362096], [146.663327264593676, -43.580853773778557], [146.04837772032036, -43.549744561538887], [145.431929559510621, -42.693776137056275], [145.295090366801702, -42.033609714527572], [144.718071323830685, -41.162551771815799], [144.743754510679679, -40.703975111657712], [145.397978143494896, -40.79254851660599]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[173.020374790740817, -40.919052422856424], [173.24723432850223, -41.331998793300784], [173.958405389702847, -40.926700534835632], [174.247586704808299, -41.349155368821755], [174.248516880589449, -41.770008233406756], [173.876446568088085, -42.233184096038826], [173.222739699595849, -42.970038344088664], [172.711246372770773, -43.372287693048591], [173.080112746470235, -43.853343601253599], [172.308583612352521, -43.865694268571346], [171.45292524646365, -44.242518812843741], [171.185137974327262, -44.897104180684892], [170.616697219116617, -45.908928724959708], [169.831422154009346, -46.355774834987614], [169.332331170934282, -46.641235446967855], [168.411353794628639, -46.619944756863603], [167.763744745146852, -46.290197442409209], [166.676886021184231, -46.219917494492265], [166.509144321964669, -45.852704766626218], [167.046424188503266, -45.110941257508671], [168.303763462596891, -44.123973077166141], [168.949408807651565, -43.935819187191427], [169.667814569373178, -43.555325616226355], [170.52491987536618, -43.03168832781283], [171.125089960004033, -42.512753594737873], [171.569713983443222, -41.767424411792135], [171.948708937871942, -41.514416599291167], [172.097227004278778, -40.956104424809766], [172.798579543344005, -40.493962090823473], [173.020374790740817, -40.919052422856424]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[174.612008905330555, -36.156397393540558], [175.336615838927202, -37.20909799575827], [175.357596470437699, -36.526193943021212], [175.808886753642668, -36.798942152657688], [175.958490025127702, -37.555381768546155], [176.763195428776754, -37.881253350578703], [177.438813104560694, -37.961248467766495], [178.010354445708884, -37.579824721020216], [178.517093540762829, -37.695373223624799], [178.27473107331403, -38.582812595373198], [177.97046023997936, -39.166342868812976], [177.206992629299322, -39.145775648760846], [176.939980503647206, -39.449736423501577], [177.032946405340311, -39.879942722331492], [176.885823602605257, -40.065977878582174], [176.508017206119376, -40.604808038089587], [176.012440220440311, -41.289624118821507], [175.239567499082995, -41.688307793953257], [175.067898391009436, -41.425894870775167], [174.650972935278645, -41.281820977545451], [175.227630243223672, -40.459235528323404], [174.900156691790158, -39.90893320084723], [173.824046665744191, -39.508854262043513], [173.852261997775514, -39.146602471677468], [174.574801874080407, -38.797683200842755], [174.743473749081062, -38.02780771255847], [174.697016636450797, -37.381128838857961], [174.292028436579386, -36.711092217761546], [174.319003534235577, -36.534823907213905], [173.840996535535822, -36.121980889634116], [173.054171177459608, -35.237125339500437], [172.636005487353799, -34.529106540669474], [173.007042271209485, -34.450661716450355], [173.551298456107673, -35.00618336358805], [174.329390497126241, -35.265495700828623], [174.612008905330555, -36.156397393540558]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[167.120011428086883, -22.159990736583524], [166.740034621444892, -22.399976088146985], [166.189732293968717, -22.129708347260475], [165.474375441752187, -21.679606621998289], [164.82981530177571, -21.149819838141951], [164.167995233413734, -20.444746595951642], [164.029605747736099, -20.105645847252347], [164.45996707586275, -20.120011895429542], [165.020036249042136, -20.459991143477737], [165.460009393575177, -20.800022067958238], [165.77998986232646, -21.080004978115596], [166.599991489933842, -21.700018812753541], [167.120011428086883, -22.159990736583524]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[178.373600000000209, -17.33992], [178.718060000000179, -17.62846], [178.55271, -18.15059], [177.932660000000197, -18.28799], [177.38146, -18.16432], [177.28504000000018, -17.72465], [177.670870000000122, -17.38114], [178.125570000000181, -17.50481], [178.373600000000209, -17.33992]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[179.36414266196428, -16.801354076946851], [178.725059362997257, -17.012041674368035], [178.596838595117248, -16.63915], [179.096609362997327, -16.433984277547438], [179.413509362997303, -16.379054277547397], [180.000000000000142, -16.06713266364244], [180.000000000000142, -16.55521656663916], [179.36414266196428, -16.801354076946851]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-179.917369384765294, -16.501783135649362], [-180.0, -16.55521656663916], [-180.0, -16.06713266364244], [-179.793320109048608, -16.020882256741231], [-179.917369384765294, -16.501783135649362]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[167.844876743845106, -16.466333103097156], [167.515181105822904, -16.59784962327997], [167.18000776597782, -16.15999521247096], [167.216801385769628, -15.891846205308454], [167.844876743845106, -16.466333103097156]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[167.107712437201513, -14.933920179913954], [167.270028111030257, -15.740020847234874], [167.001207310247906, -15.614602146062495], [166.793157993840936, -15.668810723536723], [166.649859247095577, -15.392703545801197], [166.629136997746485, -14.626497084209603], [167.107712437201513, -14.933920179913954]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[50.05651085795725, -13.555761407122006], [50.217431268114154, -14.758788750876846], [50.476536899625586, -15.226512139550593], [50.377111443895984, -15.706069431219134], [50.20027469259324, -16.000263360256824], [49.860605503138714, -15.414252618066911], [49.672606642460977, -15.7102035458025], [49.863344354050213, -16.451036879138755], [49.774564243372737, -16.875042006093651], [49.498612094934089, -17.106035658438316], [49.435618523970305, -17.953064060134423], [49.041792433474001, -19.118781019774474], [48.548540887248151, -20.49688811613413], [47.930749139198696, -22.391501153251085], [47.547723423051394, -23.781958916928488], [47.095761346226681, -24.941629733990482], [46.2824776548172, -25.178462823184148], [45.40950768411048, -25.601434421493082], [44.833573846217604, -25.346101169538926], [44.039720493349733, -24.988345228782279], [43.763768344911199, -24.460677178649973], [43.697777540874512, -23.574116306250573], [43.345654331237682, -22.776903985283866], [43.254187046081057, -22.057413018484169], [43.433297560404782, -21.336475111580185], [43.89368289569299, -21.163307386970146], [43.896370070172139, -20.830459486578164], [44.374325392439715, -20.072366224856353], [44.464397413924502, -19.435454196859055], [44.232421909366252, -18.961994724200878], [44.042976108584213, -18.331387220943185], [43.963084344261034, -17.409944756746754], [44.312468702986422, -16.850495700754919], [44.446517368351493, -16.216219170804536], [44.944936557806585, -16.179373874580435], [45.502731967964991, -15.974373467678589], [45.872993605336291, -15.793454278224672], [46.312243279817238, -15.780018405828855], [46.882182651564364, -15.210182386946343], [47.705129835812414, -14.594302666891764], [48.005214878131284, -14.091232598530382], [47.869047479042223, -13.663868503476635], [48.293827752481462, -13.784067884987479], [48.845060255738844, -13.089174899958692], [48.863508742067125, -12.487867933810477], [49.194651320193401, -12.040556735891968], [49.543518914595808, -12.469832858940592], [49.808980747279207, -12.895284925999562], [50.05651085795725, -13.555761407122006]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[143.56181115130002, -13.763655694232213], [143.922099237238939, -14.548310642152003], [144.563713820574861, -14.171176039285882], [144.894908075133486, -14.594457696188641], [145.374723748963476, -14.984976495018373], [145.271991001567329, -15.428205254785794], [145.485259637635778, -16.285672295804773], [145.637033319276981, -16.784918308176614], [145.8889042502677, -16.90692636481765], [146.160308872664558, -17.761654554925244], [146.06367394427869, -18.280072523677319], [146.387478469019669, -18.958274021075908], [147.471081577747924, -19.480722751546779], [148.177601760042506, -19.955939222902771], [148.848413527623251, -20.391209812097259], [148.717465448195611, -20.633468926681616], [149.289420200802141, -21.260510756111103], [149.678337030230665, -22.342511895438392], [150.077382440388618, -22.122783705333319], [150.482939081015189, -22.556142266533016], [150.727265252891215, -22.402404880464658], [150.899554478152339, -23.462236830338682], [151.60917524638424, -24.076256198830762], [152.073539666959135, -24.457886651306197], [152.855197381805937, -25.267501316023029], [153.136162144176836, -26.07117319102619], [153.161948683890415, -26.641319268502443], [153.092908970348589, -27.260299574494525], [153.569469028944155, -28.110066827102116], [153.512108189100246, -28.995077406532758], [153.339095493787084, -29.458201592732451], [153.069241164358885, -30.350240166954819], [153.089601678681817, -30.923641859665452], [152.891577590139406, -31.640445651986056], [152.450002476205412, -32.550002536755244], [151.709117466436766, -33.04134205498643], [151.343971795862416, -33.816023451473853], [151.010555454715217, -34.310360202777986], [150.714139439089053, -35.17345997491681], [150.328219842733262, -35.67187916437193], [150.075212030232279, -36.420205580390515], [149.946124302367224, -37.109052422841245], [149.997283970336156, -37.425260512035237], [149.423882277625495, -37.772681166333484], [148.304622430615922, -37.809061374666982], [147.381733026315288, -38.21921721776755], [146.922122837511353, -38.606532077795123], [146.317921991154861, -39.035756524411454], [145.489652134380577, -38.593767999019065], [144.876976353128185, -38.417448012039216], [145.032212355732923, -37.896187839510986], [144.485682407814039, -38.085323581699271], [143.609973586196105, -38.809465427405328], [142.745426873952994, -38.538267510737526], [142.178329705981952, -38.380034275059856], [141.606581659104762, -38.308514092767879], [140.638578729413297, -38.019332777662569], [139.992158237874349, -37.402936293285109], [139.806588169514072, -36.643602797188279], [139.574147577065304, -36.138362318670673], [139.082808058834104, -35.732754001611781], [138.120747918856324, -35.612296237939404], [138.449461704665026, -35.127261244447894], [138.207564325106688, -34.384722588845932], [137.719170363516213, -35.076825046531027], [136.829405552314739, -35.260534763328621], [137.352371047108562, -34.7073385556441], [137.503886346588359, -34.130267836240776], [137.890116001537677, -33.640478610978434], [137.81032759007914, -32.900007012668112], [136.996837192940433, -33.752771498348636], [136.37206912653167, -34.094766127256278], [135.989043410384369, -34.890118096660501], [135.208212518454133, -34.478670342752608], [135.239218377829189, -33.947953383115077], [134.613416782774635, -33.222778008763143], [134.085903761939193, -32.84807219821478], [134.273902622617101, -32.617233575166978], [132.990776808809841, -32.011224053680195], [132.288080682504955, -31.982646986622768], [131.326330601120929, -31.495803318001066], [129.535793898639696, -31.590422865527486], [128.240937534702283, -31.948488864877874], [127.102867466338324, -32.282266941051049], [126.148713820501229, -32.215966078420607], [125.088623488465686, -32.728751316052836], [124.221647983904944, -32.959486586236068], [124.028946567888596, -33.483847344701715], [123.659666782730739, -33.890179131812744], [122.811036411633694, -33.914467054989927], [122.183064406422801, -34.003402194964224], [121.299190708502607, -33.821036065406233], [120.58026818245807, -33.93017669040664], [119.89369510302825, -33.976065362281815], [119.29889936734881, -34.509366143533967], [119.00734093635802, -34.46414926527855], [118.505717808100798, -34.7468193499151], [118.024971958489544, -35.064732761374714], [117.295507440257467, -35.025458672832869], [116.625109084135005, -35.02509693780685], [115.564346958479717, -34.386427911111554], [115.026808709779601, -34.19651702243894], [115.048616164206834, -33.623425388322033], [115.545123325667163, -33.487257989232958], [115.714673700016732, -33.259571628554951], [115.679378696761404, -32.900368747694131], [115.801645135563973, -32.205062351207033], [115.689610630355247, -31.612437025683789], [115.160909051577022, -30.601594333622458], [114.997043084779506, -30.030724786094169], [115.040037876446348, -29.461095472940798], [114.641974318502065, -28.810230808224716], [114.616497837382184, -28.516398614213045], [114.173579136208474, -28.118076674107328], [114.048883905088218, -27.334765313427141], [113.477497593236961, -26.543134047147902], [113.338953078262506, -26.116545098578484], [113.778357782040274, -26.549025160429181], [113.440962355606615, -25.621278171493159], [113.936901076311671, -25.911234633082898], [114.23285200404726, -26.298446140245886], [114.216160516416977, -25.786281019801123], [113.72125532435777, -24.998938897402127], [113.625343866024053, -24.683971042583153], [113.393523390762681, -24.384764499613269], [113.502043898575579, -23.806350192970257], [113.706992629045175, -23.560215345964068], [113.843418410295698, -23.059987481378741], [113.736551548316157, -22.475475355725379], [114.149756300922007, -21.755881036061012], [114.225307244932679, -22.517488295178723], [114.647762078918817, -21.829519952077007], [115.460167270979326, -21.495173435148558], [115.947372674627019, -21.068687839443712], [116.711615431791557, -20.701681817306834], [117.166316359527713, -20.623598728113819], [117.441545037914324, -20.746898695562251], [118.22955895393298, -20.374208265873236], [118.836085239742744, -20.263310642174844], [118.987807244951767, -20.044202569257322], [119.252493931150724, -19.952941989829839], [119.805225050944586, -19.976506442954985], [120.856220330896718, -19.683707777589191], [121.399856398607227, -19.239755547769732], [121.655137974129076, -18.705317885007133], [122.241665480641842, -18.197648614171854], [122.28662397673574, -17.798603204014015], [122.312772251475423, -17.254967136303463], [123.012574497571933, -16.405199883695872], [123.433789097183109, -17.268558037996229], [123.859344517106678, -17.069035332917267], [123.503242222183275, -16.596506036040466], [123.817073195491872, -16.111316013251994], [124.258286574399875, -16.327943617419564], [124.379726190285822, -15.567059828353976], [124.926152785340122, -15.075100192935409], [125.16727501841396, -14.680395603090091], [125.670086704613851, -14.510070082256021], [125.685796340030521, -14.230655612853838], [126.12514936737611, -14.347340996968953], [126.142822707219892, -14.095986830301229], [126.582589146023764, -13.952791436420497], [127.065867140817346, -13.817967624570926], [127.804633416861947, -14.276906019755131], [128.359689976108967, -14.869169610252271], [128.985543247595928, -14.875990899314758], [129.621473423379626, -14.969783623924556], [129.409600050983016, -14.420669854391122], [129.888640578328619, -13.618703301653497], [130.339465773642956, -13.357375583553477], [130.183506300986068, -13.107520033422304], [130.617795037967056, -12.536392103732481], [131.223494500860028, -12.183648776908214], [131.735091180549517, -12.302452894747162], [132.575298293183124, -12.114040622611014], [132.55721154188106, -11.603012383676699], [131.824698114143729, -11.2737818335452], [132.357223748911423, -11.128519382372744], [133.019560581596437, -11.376411228076847], [133.550845981989113, -11.786515394745138], [134.393068475482011, -12.042365411022175], [134.678632440326993, -11.941182956594702], [135.298491245668032, -12.248606052299053], [135.882693312727696, -11.962266940969798], [136.258380975489473, -12.049341729381609], [136.492475213771655, -11.857208754120393], [136.951620314685016, -12.351958916882836], [136.685124953355768, -12.887223402562057], [136.305406528875181, -13.291229750219898], [135.961758254134196, -13.324509372615893], [136.077616815332561, -13.724278252825783], [135.783836297753254, -14.223989353088214], [135.428664178611228, -14.7154322241839], [135.500184360903262, -14.997740573794445], [136.295174595281395, -15.550264987859123], [137.065360142159506, -15.870762220933358], [137.580470819244823, -16.215082289294088], [138.303217401279056, -16.807604261952761], [138.585164015863398, -16.806622409739177], [139.108542922115561, -17.062679131745369], [139.260574985918225, -17.371600843986187], [140.215245396078245, -17.710804945550066], [140.87546349503927, -17.369068698803943], [141.071110467696286, -16.832047214426723], [141.274095493738827, -16.388870131091693], [141.398222284103866, -15.840531508042588], [141.702183058844668, -15.044921156476931], [141.563380161708693, -14.561333103089609], [141.635520461188122, -14.270394789286286], [141.519868605718983, -13.698078301653823], [141.650920038011094, -12.944687595270565], [141.842691278246292, -12.74154753993129], [141.686990187750865, -12.407614434461152], [141.928629185147571, -11.877465915578796], [142.118488397388063, -11.328042087451621], [142.143706496346368, -11.042736504768243], [142.515260044524979, -10.668185723516729], [142.797310011974076, -11.157354831591533], [142.866763136974299, -11.784706719614931], [143.115946893485699, -11.905629571177911], [143.158631626558787, -12.325655612846205], [143.522123651299893, -12.834358412327433], [143.597157830987641, -13.400422051652598], [143.56181115130002, -13.763655694232213]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[162.11902469304087, -10.482719008021135], [162.398645868172224, -10.826367282762121], [161.700032180018383, -10.820011081590238], [161.319796991214815, -10.204751478723225], [161.917383254238047, -10.446700534713756], [162.11902469304087, -10.482719008021135]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[120.715608758630452, -10.239581394087878], [120.295014276206899, -10.258649997603627], [118.967808465654713, -9.557969252158045], [119.900309686361567, -9.361340427287516], [120.425755649905369, -9.665921319215798], [120.775501743656804, -9.969675388227458], [120.715608758630452, -10.239581394087878]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[160.85222863183796, -9.872937106977105], [160.462588332357285, -9.895209649294841], [159.849447463214204, -9.794027194867368], [159.640002883135224, -9.639979750205271], [159.702944777666715, -9.242949720906779], [160.362956170898457, -9.400304457235549], [160.688517694337264, -9.610162448772911], [160.85222863183796, -9.872937106977105]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[161.679981724289149, -9.599982191611375], [161.529396600590559, -9.784312025596449], [160.788253208660564, -8.91754322676492], [160.579997186524423, -8.320008640173967], [160.920028111004939, -8.320008640173967], [161.280006138349989, -9.120011488484451], [161.679981724289149, -9.599982191611375]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[124.435950148619412, -10.140000909061442], [123.579981724136672, -10.359987481327963], [123.459989048355027, -10.239994805546189], [123.550009393407464, -9.90001555749798], [123.980008986508125, -9.290026950724695], [124.968682489116276, -8.892790215697133], [125.086246372580291, -8.656887302284673], [125.947072381698348, -8.432094821815014], [126.64470421763852, -8.398246758663859], [126.957243280139807, -8.273344821814362], [127.335928175974601, -8.397316582882652], [126.967991978056517, -8.668256117388935], [125.925885044458681, -9.10600717533336], [125.088520135601101, -9.393173109579322], [124.435950148619412, -10.140000909061442]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[117.900018345207769, -8.095681247594925], [118.260616489740499, -8.362383314653329], [118.878459914222134, -8.28068287519983], [119.126506789223157, -8.705824883665073], [117.970401645989313, -8.906639499551346], [117.277730747549043, -9.040894870645573], [116.740140822416691, -9.032936700072639], [117.083737420725328, -8.457157891476626], [117.632024367342154, -8.449303073768192], [117.900018345207769, -8.095681247594925]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[122.903537225436111, -8.094234307490751], [122.756982863456329, -8.649807631060739], [121.254490594570171, -8.933666273639943], [119.924390903809609, -8.810417982623875], [119.920928582846074, -8.444858900591171], [120.715091994307642, -8.236964613480964], [121.341668735846582, -8.536739597206108], [122.007364536630433, -8.460620212440162], [122.903537225436111, -8.094234307490751]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[159.875027297198613, -8.337320244991716], [159.917401971677947, -8.538289890174866], [159.133677199539392, -8.1141814103554], [158.586113722974773, -7.754823500197716], [158.211149530264862, -7.421872246941248], [158.359977655265453, -7.320017998893931], [158.820001255527785, -7.560003350457393], [159.640002883135224, -8.020026950719668], [159.875027297198613, -8.337320244991716]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[157.538425734689241, -7.347819919466943], [157.339419793933274, -7.404767347852655], [156.902030471014854, -7.176874281445407], [156.491357863591332, -6.765943291860495], [156.542827590154019, -6.59933847415148], [157.140000441718911, -7.021638278840655], [157.538425734689241, -7.347819919466943]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[108.623478631629013, -6.777673841990691], [110.539227329553313, -6.877357679881698], [110.759575636845938, -6.465186455921753], [112.614811232556434, -6.946035658397605], [112.978768345188115, -7.59421314863458], [114.478935174621171, -7.776527601760293], [115.705526971501143, -8.370806573116866], [114.564511346496516, -8.751816908404834], [113.464733514460903, -8.348947442257426], [112.559672479301042, -8.376180922075264], [111.522061395312534, -8.302128594600958], [110.586149530074323, -8.122604668819022], [109.427667270955197, -7.740664157749762], [108.693655226681386, -7.641600437046222], [108.277763299596387, -7.766657403192582], [106.454102004016164, -7.354899590690962], [106.280624220812371, -6.924899997590302], [105.365486281355544, -6.851416110871256], [106.051645949327082, -5.8959188777945], [107.265008579540194, -5.954985039904059], [108.072091099074697, -6.345762220895253], [108.486846144649377, -6.421984958525769], [108.623478631629013, -6.777673841990691]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[134.724624465066682, -6.214400730009302], [134.210133905168931, -6.895237725454706], [134.112775506731026, -6.142467136259015], [134.290335728085864, -5.783057549669039], [134.49962527886791, -5.445042006047899], [134.727001580952134, -5.73758228925216], [134.724624465066682, -6.214400730009302]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[155.88002566957843, -6.81999684003776], [155.599991082988794, -6.919990736522493], [155.166994256815173, -6.535931491729315], [154.729191522438356, -5.900828138862209], [154.514114211239729, -5.139117526880014], [154.652503696917307, -5.042430922061939], [154.75999067608447, -5.339983819198494], [155.062917922179366, -5.566791680527501], [155.547746209941778, -6.200654799019659], [156.019965448224838, -6.540013929880402], [155.88002566957843, -6.81999684003776]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[151.982795851854547, -5.478063246282346], [151.459106887008687, -5.560280450058741], [151.301390415653856, -5.840728448106788], [150.754447056276689, -6.083762709175474], [150.241196730753899, -6.317753594593086], [149.709963006793345, -6.316513360218053], [148.89006473205049, -6.026040134305433], [148.318936802360753, -5.747142429226145], [148.401825799756892, -5.437755629094724], [149.298411900020852, -5.583741550319317], [149.845561965127303, -5.505503431829339], [149.996250441690307, -5.026101169457675], [150.139755894164892, -5.001348158389888], [150.236907586873571, -5.532220147324281], [150.807467075808148, -5.455842380396888], [151.089672072554009, -5.113692722192383], [151.647880894170925, -4.757073662946183], [151.537861769821546, -4.167807305521976], [152.136791620084381, -4.14879037843852], [152.338743117481016, -4.312966403829861], [152.318692661751783, -4.867661228050764], [151.982795851854547, -5.478063246282346]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[127.249215122588993, -3.45906503663889], [126.874922723498884, -3.79098276124958], [126.18380211802733, -3.607376397316571], [125.989033644719285, -3.177273451351326], [127.000651483265045, -3.129317722184496], [127.249215122588993, -3.45906503663889]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[130.471344028851803, -3.09376433676762], [130.834836053592852, -3.858472181822762], [129.990546502808144, -3.446300957862832], [129.155248651242431, -3.362636813982263], [128.590683628453661, -3.428679294451257], [127.898891229362363, -3.393435967628207], [128.135879347852807, -2.843650404475014], [129.370997756060973, -2.802154229344652], [130.471344028851803, -3.09376433676762]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[153.140037876598768, -4.499983412294114], [152.827292108368312, -4.766427097190999], [152.638673130503008, -4.176127211120928], [152.406025832324957, -3.789742526874562], [151.953236932583621, -3.462062269711822], [151.384279413050052, -3.035421644710112], [150.662049595338857, -2.741486097833956], [150.939965448204561, -2.500002129734028], [151.479984165654599, -2.779985039891386], [151.820015090135115, -2.999971612157907], [152.239989455371102, -3.240008640153661], [152.640016717742554, -3.659983005389748], [153.019993524384716, -3.980015150573294], [153.140037876598768, -4.499983412294114]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[134.1433679546478, -1.151867364103595], [134.422627394753107, -2.769184665542397], [135.457602980694759, -3.367752780779128], [136.293314243718811, -2.30704233155619], [137.440737746327528, -1.703513278819372], [138.329727411044786, -1.70268645590275], [139.184920689042968, -2.051295668143638], [139.926684198160416, -2.409051608900285], [141.000210402591875, -2.600151055515639], [142.7352466167915, -3.289152927263217], [144.583970982033264, -3.861417738463402], [145.273179559509998, -4.373737888205028], [145.829786411725735, -4.876497897972683], [145.981921828392984, -5.465609226100028], [147.648073358347602, -6.083659356310818], [147.891107619416204, -6.61401458092233], [146.97090538959489, -6.721656589386356], [147.191873814074967, -7.388024183789994], [148.084635858349401, -8.044108168167611], [148.734105259393601, -9.104663588093771], [149.306835158484461, -9.071435642130069], [149.266630894161352, -9.514406019736029], [150.038728469034282, -9.684318129111702], [149.73879845601229, -9.872937106977105], [150.801627638959161, -10.293686618697521], [150.690574985963877, -10.582712904505968], [150.028393182575797, -10.652476088099945], [149.782310012002, -10.393267103723943], [148.923137648717301, -10.280922539921363], [147.913018426708021, -10.130440769087471], [147.135443150012264, -9.492443536012019], [146.567880894150591, -8.942554619994155], [146.048481073184945, -8.06741423913131], [144.744167922138075, -7.630128269077474], [143.89708784400969, -7.915330498896282], [143.286375767184353, -8.245491224809072], [143.413913202080693, -8.983068942911032], [142.628431431244252, -9.326820570516503], [142.068258905200281, -9.159595635620036], [141.033851760013903, -9.117892754760518], [140.143415155192628, -8.297167657100957], [139.127766554928115, -8.096042982621029], [138.881476678625035, -8.38093515384611], [137.614473911692841, -8.411682631059762], [138.03909915583526, -7.597882175327356], [138.668621454014811, -7.320224704623087], [138.407913853102372, -6.232849216337485], [137.927839797110863, -5.393365573756], [135.989250116113482, -4.546543877789063], [135.164597609599724, -4.462931410340872], [133.662880487197953, -3.538853448097527], [133.367704705946807, -4.024818617370315], [132.983955519747354, -4.112978610860281], [132.756940952689007, -3.74628264731713], [132.753788690319283, -3.311787204607072], [131.989804315316206, -2.820551039240556], [133.066844517143494, -2.460417982598443], [133.780030959203572, -2.47984832114021], [133.696211786026169, -2.214541517753688], [132.232373488494289, -2.21252613689434], [131.836221958544769, -1.617161960459697], [130.942839797082883, -1.432522067880811], [130.519558140180067, -0.937720228686075], [131.867537876513637, -0.695461114101818], [132.380116408416797, -0.369537855636977], [133.985548130428441, -0.780210463060456], [134.1433679546478, -1.151867364103595]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[125.240500522971587, 1.419836127117605], [124.437035353697439, 0.427881171058957], [123.68550499887678, 0.235593166500877], [122.723083123872897, 0.431136786293337], [121.056724888189109, 0.381217352699352], [120.183083123862758, 0.237246812334206], [120.040869582195484, -0.519657891444865], [120.935905389490785, -1.408905938323372], [121.475820754076182, -0.955962009285116], [123.340564813328484, -0.615672702643181], [123.258399285984495, -1.076213067228338], [122.822715285331668, -0.930950616055881], [122.388529901215378, -1.516858005381124], [121.50827355355554, -1.904482924002423], [122.454572381684358, -3.186058444840967], [122.271896193532569, -3.529500013852712], [123.170962762546566, -4.683693129091722], [123.162332798353788, -5.340603936386046], [122.628515252778726, -5.634591159694494], [122.236394484548072, -5.282933037948297], [122.719569126477069, -4.46417164471589], [121.738233677254442, -4.851331475446585], [121.489463332201325, -4.574552504091315], [121.619171177253946, -4.188477878438675], [120.898181593917712, -3.602105401222829], [120.972388950688782, -2.62764291749491], [120.305452915529912, -2.931603692235726], [120.390047235191759, -4.097579034037309], [120.430716587405385, -5.528241062037779], [119.796543410319572, -5.673400160345651], [119.366905552244958, -5.379878024927805], [119.653606398600203, -4.459417412944958], [119.498835483885983, -3.494411716326525], [119.078344354327015, -3.487021986508765], [118.76776899625284, -2.801999200047689], [119.180973748858747, -2.147103773612798], [119.323393996255135, -1.353147067880485], [119.825998976725913, 0.154254462073482], [120.035701938966355, 0.566477362465704], [120.885779250167644, 1.309222723796836], [121.666816847826993, 1.013943589681077], [122.927566766451861, 0.875192368977366], [124.077522414242907, 0.917101955566139], [125.065989211121831, 1.643259182131544], [125.240500522971587, 1.419836127117605]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[128.688248732620735, 1.132385972494006], [128.63595218314137, 0.258485826006179], [128.120169712436194, 0.356412665199272], [127.968034295768888, -0.252077325037533], [128.37999881399972, -0.780003757331301], [128.100015903842376, -0.899996433113074], [127.696474644075039, -0.266598402511505], [127.399490187693772, 1.011721503092559], [127.600511509309086, 1.810690822757181], [127.93237755748757, 2.174596258956555], [128.004156121940838, 1.628531398928331], [128.594559360875479, 1.540810655112864], [128.688248732620735, 1.132385972494006]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[105.81765506390937, -5.852355645372413], [104.710384149191526, -5.873284600450646], [103.868213332130807, -5.037314955264989], [102.584260695406982, -4.220258884298204], [102.156173130301084, -3.61414600994685], [101.399113397225079, -2.799777113459186], [100.902502882900222, -2.05026213949786], [100.141980828860682, -0.650347588710972], [99.263739862060305, 0.183141587724649], [98.970011020913347, 1.042882391764536], [98.601351352943112, 1.823506577965532], [97.69959760944991, 2.453183905442017], [97.176942173249898, 3.30879059489861], [96.424016554757287, 3.868859768077911], [95.380876092513546, 4.97078217205366], [95.293026157617334, 5.479820868344817], [95.936862827541773, 5.439513251157095], [97.484882033277103, 5.246320909033912], [98.369169142655693, 4.268370266126368], [99.142558628335934, 3.59034963624083], [99.693997837322428, 3.174328518075157], [100.641433546961679, 2.099381211755698], [101.658012323007455, 2.083697414555189], [102.49827111207324, 1.398700466310203], [103.076840448013144, 0.56136139566884], [103.838396030698362, 0.104541734208667], [103.437645298274987, -0.711945896002945], [104.010788608824072, -1.059211521004329], [104.369991489684963, -1.084843031421116], [104.53949018760224, -1.782371514496802], [104.887892694114072, -2.340425306816755], [105.622111444116996, -2.42884368246807], [106.108593377712708, -3.06177662517895], [105.857445916774196, -4.305524997579809], [105.81765506390937, -5.852355645372413]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[117.875627069165972, 1.827640692548897], [118.996747267738186, 0.902219143066048], [117.811858351717802, 0.784241848143722], [117.478338657706075, 0.102474676917012], [117.521643507966616, -0.80372323975331], [116.560048455879524, -1.487660821136231], [116.533796828275257, -2.483517347832901], [116.148083937648636, -4.012726332214015], [116.00085778204911, -3.657037448749108], [114.864803094544556, -4.106984144714417], [114.468651564595092, -3.495703627133835], [113.755671828264127, -3.43916961020652], [113.256994256647573, -3.118775729996955], [112.068126255340786, -3.478392022316086], [111.703290643360077, -2.994442233902646], [111.048240187628238, -3.049425957861203], [110.223846063276, -2.934032484553484], [110.070935500124364, -1.592874037282499], [109.571947869914055, -1.314906507984489], [109.091873813922547, -0.459506524257151], [108.952657505328233, 0.415375474444346], [109.069136183714107, 1.341933905437642], [109.663260125773803, 2.006466986494985], [110.396135288537124, 1.663774725751381], [111.168852980597507, 1.850636704918784], [111.370081007942105, 2.697303371588859], [111.796928338672927, 2.885896511238073], [112.995614862115275, 3.102394924324869], [113.712935418758747, 3.893509426281128], [114.204016554828428, 4.52587392823682], [114.599961379048722, 4.900011298029924], [115.450710483869813, 5.447729803891562], [116.220741001451046, 6.143191229675523], [116.725102980619766, 6.924771429873999], [117.129626092600489, 6.928052883324568], [117.643393182446317, 6.422166449403207], [117.689075148592309, 5.987490139180181], [118.347691278152269, 5.708695786965464], [119.181903924639954, 5.407835598162151], [119.110693800941789, 5.016128241389765], [118.43972700406411, 4.966518866389606], [118.618320754064854, 4.478202419447541], [117.882034946770176, 4.137551377779488], [117.313232456533541, 3.234428208830579], [118.048329705885436, 2.287690131027361], [117.875627069165972, 1.827640692548897]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[126.37681359263749, 8.414706325713254], [126.478512811387901, 7.750354112168978], [126.537423944200583, 7.189380601424574], [126.196772902532558, 6.274294338400026], [125.831420526229095, 7.293715318221842], [125.363852166852297, 6.786485297060892], [125.683160841983721, 6.049656887227243], [125.396511672060655, 5.58100332277229], [124.219787632342417, 6.161355495626168], [123.938719517106932, 6.885135606306122], [124.24366214406129, 7.360610459823661], [123.61021243702757, 7.833527329942754], [123.296071405125275, 7.418875637232773], [122.825505812675459, 7.457374579290203], [122.085499302255727, 6.899424139834835], [121.919928013192617, 7.192119452335973], [122.31235884001714, 8.034962063016408], [122.942397902519673, 8.316236883981091], [123.487687616063539, 8.693009751821194], [123.841154412939858, 8.240324204944386], [124.601469761250229, 8.514157619659017], [124.764612257995708, 8.96040945071546], [125.471390822451639, 8.986996975129628], [125.412117954612853, 9.760334784377548], [126.222714471543185, 9.286074327018838], [126.306636997585173, 8.782487494334561], [126.37681359263749, 8.414706325713254]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[81.218019647144303, 6.197141424988303], [80.348356968104468, 5.968369859232141], [79.87246870312859, 6.763463446474915], [79.695166863935157, 8.200843410673372], [80.147800734379672, 9.824077663609557], [80.838817986986641, 9.268426825391174], [81.304319289071799, 8.564206244333675], [81.787959018891428, 7.523055324733178], [81.637322218760659, 6.481775214051936], [81.218019647144303, 6.197141424988303]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-60.935, 10.11], [-61.77, 10.0], [-61.95, 10.09], [-61.66, 10.365], [-61.68, 10.76], [-61.104999999999876, 10.89], [-60.895, 10.855], [-60.935, 10.11]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[123.982437778825869, 10.278778591345713], [123.623183221532798, 9.950090643753299], [123.309920688979417, 9.318268744336677], [122.99588300994165, 9.0221886255204], [122.380054966319477, 9.713360907424203], [122.586088901867157, 9.981044826696106], [122.837081333508792, 10.261156927934238], [122.947410516451924, 10.881868394408031], [123.498849725438475, 10.940624497923949], [123.33777428598475, 10.267383938025347], [124.077935825701246, 11.23272553145371], [123.982437778825869, 10.278778591345713]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[118.504580926590364, 9.316382554558004], [117.174274530100746, 8.367499904814665], [117.664477166821399, 9.06688873945285], [118.386913690261764, 9.684499619989225], [118.987342157061079, 10.376292019080495], [119.511496209797571, 11.369668077027214], [119.689676548339975, 10.554291490109861], [119.029458449379007, 10.003653265823772], [118.504580926590364, 9.316382554558004]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[121.883547804859091, 11.89175507247198], [122.483821242361529, 11.582187404827508], [123.12021650603603, 11.58366018314787], [123.100837843926541, 11.165933742716476], [122.637713657726721, 10.741308498574128], [122.002610304859587, 10.441016750526089], [121.967366978036551, 10.905691229694625], [122.038370396005547, 11.415840969279941], [121.883547804859091, 11.89175507247198]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[125.502551711123516, 12.162694606978249], [125.78346479706218, 11.046121934447754], [125.011883986512288, 11.31145457605038], [125.0327612651582, 10.975816148314706], [125.277449172060216, 10.358722032101312], [124.801819289245742, 10.134678859899893], [124.760168084818503, 10.837995103392203], [124.459101190286077, 10.889929917845535], [124.302521600441736, 11.49537099857713], [124.891012811381557, 11.415582587118493], [124.87799035044398, 11.794189968304906], [124.266761509295662, 12.557760931849685], [125.227116327007849, 12.535720933477194], [125.502551711123516, 12.162694606978249]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[121.52739383350351, 13.069590155484519], [121.262190382981572, 12.205560207564403], [120.833896112146618, 12.704496161342419], [120.323436313967505, 13.46641347905377], [121.180128208502111, 13.429697373910443], [121.52739383350351, 13.069590155484519]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[121.321308221523594, 18.504064642810903], [121.937601353036456, 18.218552354398469], [122.246006300954292, 18.478949896717097], [122.33695682178805, 18.224882717354063], [122.174279412933259, 17.810282701076375], [122.515653924653378, 17.093504746971973], [122.252310825693911, 16.262444362854012], [121.662786086108241, 15.9310175643501], [121.505069614753467, 15.124813544164624], [121.728828566577278, 14.328376369682246], [122.258925409027341, 14.218202216035976], [122.701275669445721, 14.336541245984336], [123.950295037940322, 13.78213064214097], [123.855107049658642, 13.237771104378382], [124.181288690284902, 12.997527370653472], [124.077419061378322, 12.536676947474575], [123.298035109552274, 13.027525539598898], [122.928651971530002, 13.552919826710408], [122.671355015148691, 13.185836289925049], [122.034649692880549, 13.784481919810247], [121.126384718918615, 13.636687323455547], [120.628637323083325, 13.857655747935553], [120.679383579593917, 14.271015529838309], [120.991819289230563, 14.525392767794983], [120.693336216312701, 14.756670640517285], [120.564145135583061, 14.396279201713824], [120.070428501466409, 14.970869452367197], [119.920928582846074, 15.40634674729074], [119.883773228028275, 16.363704331929966], [120.28648766487882, 16.034628811095331], [120.390047235191759, 17.599081122299509], [120.715867140791971, 18.505227362537454], [121.321308221523594, 18.504064642810903]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-65.591003790942892, 18.228034979723873], [-65.847163865813656, 17.975905666571805], [-66.599934455009446, 17.981822618069288], [-67.184162360285171, 17.946553453030134], [-67.242427537694397, 18.374460150622866], [-67.100679083917754, 18.520601101144422], [-66.282434455008143, 18.514761664295321], [-65.771302863209343, 18.426679185453935], [-65.591003790942892, 18.228034979723873]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-76.902561408175615, 17.868237819891675], [-77.206341315403478, 17.701116237859793], [-77.766022915340557, 17.861597398342241], [-78.337719285785482, 18.225967922432318], [-78.217726610003893, 18.454532782459324], [-77.797364671525685, 18.524218451404721], [-77.569600796199097, 18.49052541755043], [-76.896618618462043, 18.400866807524096], [-76.365359056285428, 18.160700588447554], [-76.19965857614153, 17.886867173732924], [-76.902561408175615, 17.868237819891675]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-72.579672817663635, 19.871500555902344], [-71.712361416292936, 19.714455878167314], [-71.587304450146576, 19.884910590082058], [-70.806706102161684, 19.880285549391971], [-70.214364997016077, 19.622885240146104], [-69.950815192327553, 19.647999986239995], [-69.769250047470081, 19.293267116772483], [-69.222125820579805, 19.313214219637047], [-69.254346076113734, 19.015196234609988], [-68.809411994080705, 18.979074408437839], [-68.317943284768916, 18.612197577381636], [-68.689315965434531, 18.205142320218556], [-69.164945848248834, 18.422648423735126], [-69.623987596297525, 18.380712998930363], [-69.952933926051486, 18.428306993071061], [-70.133232998317936, 18.245915025296881], [-70.517137213814266, 18.184290879788904], [-70.669298468697576, 18.426885891182991], [-70.999950120717159, 18.283328762276341], [-71.400209927033842, 17.598564357976528], [-71.657661912711916, 17.757572740138727], [-71.708304816357952, 18.044997056546208], [-72.372476162389262, 18.214960842354088], [-72.844411180294884, 18.145611070218337], [-73.454554816365032, 18.217906398994813], [-73.922433234335557, 18.030992743394933], [-74.458033616824707, 18.342549953682663], [-74.369925299767118, 18.664907538319397], [-73.449542202432752, 18.526052964751102], [-72.694937099890666, 18.445799465401791], [-72.334881557897035, 18.668421535715311], [-72.79164954292483, 19.10162506761796], [-72.784104783810307, 19.483591416903408], [-73.41502234566164, 19.639550889560297], [-73.18979061551758, 19.915683905511997], [-72.579672817663635, 19.871500555902344]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[110.339187860151497, 18.678395087147578], [109.475209588663716, 18.197700913968617], [108.655207961056163, 18.507681993071401], [108.626217482540511, 19.367887885001863], [109.119055617308078, 19.82103851976936], [110.211598748822865, 20.101253973871991], [110.786550734502242, 20.077534491450052], [111.010051304164648, 19.695929877190821], [110.570646600386766, 19.255879218009312], [110.339187860151497, 18.678395087147578]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-155.54211, 19.08348], [-155.68817, 18.916190000000114], [-155.93665, 19.05939], [-155.90806, 19.338880000000103], [-156.07347, 19.70294], [-156.02368, 19.81422], [-155.85008, 19.97729], [-155.91907, 20.17395], [-155.86108, 20.26721], [-155.78505, 20.2487], [-155.40214, 20.07975], [-155.22452, 19.99302], [-155.06226, 19.8591], [-154.80741, 19.508710000000121], [-154.83147, 19.45328], [-155.22217, 19.239720000000119], [-155.54211, 19.08348]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-156.07926, 20.64397], [-156.41445, 20.57241], [-156.58673, 20.783], [-156.70167, 20.8643], [-156.71055, 20.92676], [-156.61258, 21.01249], [-156.25711, 20.91745], [-155.99566, 20.76404], [-156.07926, 20.64397]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-156.75824, 21.17684], [-156.78933, 21.06873], [-157.32521, 21.09777], [-157.25027, 21.21958], [-156.75824, 21.17684]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-157.65283, 21.32217], [-157.70703, 21.26442], [-157.7786, 21.27729], [-158.12667, 21.31244], [-158.2538, 21.53919], [-158.29265, 21.57912], [-158.0252, 21.71696], [-157.94161, 21.65272], [-157.65283, 21.32217]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-159.34512, 21.982], [-159.46372, 21.88299], [-159.80051, 22.06533], [-159.74877, 22.1382], [-159.5962, 22.23618], [-159.36569, 22.21494], [-159.34512, 21.982]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-79.679523688460307, 22.765303249598816], [-79.281485968732028, 22.399201565027028], [-78.347434455056401, 22.512166246017074], [-77.993295864560167, 22.277193508385892], [-77.146422492161008, 21.657851467367806], [-76.523824835908442, 21.206819566324327], [-76.194620123993133, 21.22056549731397], [-75.598222418912627, 21.016624457274077], [-75.671060350228117, 20.73509125414796], [-74.93389604358444, 20.693905137611324], [-74.178024868451246, 20.284627793859755], [-74.296648118777142, 20.050378526280696], [-74.961594611292981, 19.923435370355691], [-75.63468014189462, 19.873774318923154], [-76.323656175426009, 19.952890936762117], [-77.755480923153101, 19.855480861891891], [-77.085108405246672, 20.413353786698792], [-77.49265458851653, 20.673105373613851], [-78.137292243141502, 20.739948838783477], [-78.482826707661161, 21.028613389565805], [-78.719866502583983, 21.598113511638417], [-79.284999966127856, 21.559175319906473], [-80.217475348618535, 21.827324327068965], [-80.517534552721429, 22.037078965741756], [-81.82094336620321, 22.19205658618506], [-82.169991828118697, 22.387109279870742], [-81.795001797192583, 22.636964830002086], [-82.775897996740781, 22.688150336187107], [-83.494458787759299, 22.168517971276088], [-83.90880042187564, 22.154565334557304], [-84.052150845053319, 21.910575059491322], [-84.547030198896437, 21.801227728761575], [-84.974911058273165, 21.896028143801061], [-84.447062140627793, 22.204949856041878], [-84.230357021811841, 22.56575470630375], [-83.778239915690108, 22.788118394455637], [-83.267547573565622, 22.98304189706073], [-82.510436164057495, 23.078746649665135], [-82.268151211257049, 23.188610744717664], [-81.404457160146848, 23.117271429938768], [-80.61876868358118, 23.105980129482958], [-79.679523688460307, 22.765303249598816]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-77.534659999999889, 23.75975], [-77.78, 23.71], [-78.03405, 24.28615], [-78.40848, 24.57564], [-78.19087, 25.2103], [-77.89, 25.17], [-77.54, 24.340000000000117], [-77.534659999999889, 23.75975]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[121.175632358892727, 22.790857245367135], [120.747079705896255, 21.970571397382088], [120.220083449383679, 22.814860948166682], [120.106188592612426, 23.556262722258225], [120.694679803552305, 24.538450832613734], [121.495044386888793, 25.295458889257361], [121.951243931161542, 24.997595933526981], [121.777817824389928, 24.394273586519432], [121.175632358892727, 22.790857245367135]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-77.82, 26.58], [-78.91, 26.42], [-78.979999999999876, 26.79], [-78.51, 26.87], [-77.85, 26.84], [-77.82, 26.58]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-77.0, 26.590000000000117], [-77.17255, 25.87918], [-77.35641, 26.00735], [-77.34, 26.53], [-77.78802, 26.92516], [-77.79, 27.04], [-77.0, 26.590000000000117]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[134.638428176003885, 34.149233710256425], [134.766379022358564, 33.80633474378368], [134.203415968970916, 33.201177883429608], [133.792950067276507, 33.521985175097683], [133.280268182508877, 33.289570420864834], [133.014858026257883, 32.704567369104694], [132.363114862192759, 32.98938202568138], [132.371176385630264, 33.46364248303999], [132.924372593314814, 34.060298570282129], [133.49296837782228, 33.944620876596673], [133.904106073136376, 34.364931138642703], [134.638428176003885, 34.149233710256425]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[34.576473829900607, 35.671595567358764], [33.900804477684261, 35.245755927057587], [33.973616570783491, 35.058506374647976], [34.004880812320124, 34.978097846001873], [32.979827101378504, 34.571869411755415], [32.490296258277596, 34.701654771456532], [32.256667107885988, 35.103232326796615], [32.731780226377538, 35.140025946588423], [32.802473585752892, 35.145503648411392], [32.94696089044092, 35.386703396133697], [33.66722700372506, 35.373215847305602], [34.576473829900607, 35.671595567358764]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[23.699980096133089, 35.705004380835618], [24.246665073348765, 35.368022365860156], [25.025015496528965, 35.424995632462071], [25.76920779796427, 35.35401805270908], [25.745023227651643, 35.179997666966216], [26.290002882601726, 35.299990342748004], [26.164997592887744, 35.004995429009881], [24.72498213064236, 34.91998769788961], [24.735007358506948, 35.084990546197673], [23.51497846852817, 35.279991563450949], [23.699980096133089, 35.705004380835618]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[15.520376010813834, 38.231155096991557], [15.160242954171821, 37.44404551853782], [15.309897902089091, 37.1342194687318], [15.099988234119536, 36.619987290995397], [14.335228712632016, 36.996630967754726], [13.826732618880015, 37.104531358380115], [12.431003859108785, 37.612949937483705], [12.570943637755221, 38.126381130519604], [13.74115644700467, 38.034965521795442], [14.761249220446246, 38.143873602850505], [15.520376010813834, 38.231155096991557]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[9.210011834356351, 41.209991360024304], [9.809975213264977, 40.500008856766101], [9.669518670295645, 39.177376410471879], [9.214817742559489, 39.240473334300134], [8.806935662479702, 38.906617743478478], [8.4283024430772, 39.171847032216618], [8.388253208050941, 40.378310858718777], [8.159998406617746, 40.950007229163703], [8.709990675500194, 40.899984442705232], [9.210011834356351, 41.209991360024304]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[140.976387567305295, 37.142074286440163], [140.599769728762141, 36.343983466124456], [140.774074334882613, 35.842877102190215], [140.253279250245185, 35.138113918593746], [138.975527785396224, 34.667600002576108], [137.217598911691283, 34.606285915661772], [135.792983026268899, 33.464805202766712], [135.120982700745429, 33.84907115328906], [135.079434849182775, 34.59654490817482], [133.340316196831992, 34.375938218720847], [132.156770868051325, 33.904933376596517], [130.98614464734348, 33.885761420216198], [132.00003624891005, 33.149992377244615], [131.332790155157369, 31.450354519164847], [130.686317987186015, 31.029579169228327], [130.202419875204981, 31.418237616495389], [130.447676222862157, 32.319474595665639], [129.814691603718899, 32.610309556604363], [129.408463169472526, 33.29605581311759], [130.353935174684665, 33.604150702441785], [130.878450962447147, 34.23274282484013], [131.88422936414392, 34.749713853487918], [132.617672967662514, 35.433393052709391], [134.608300815977799, 35.73161774346579], [135.67753787652893, 35.527134100886826], [136.723830601142453, 37.304984239240383], [137.390611607004502, 36.827390651998826], [138.857602166906219, 37.827484646143546], [139.426404657142911, 38.215962225897556], [140.054790073812086, 39.43880748143647], [139.883379347899876, 40.563312486323611], [140.30578250545372, 41.195005194659643], [141.368973423426695, 41.378559882160374], [141.914263136970561, 39.991616115878685], [141.884600864835051, 39.180864569651419], [140.959489373945843, 38.17400096287659], [140.976387567305295, 37.142074286440163]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[9.560016310269219, 42.152491970379458], [9.229752231491858, 41.380006822264448], [8.775723097375447, 41.583611965494413], [8.544212680707801, 42.256516628583], [8.746009148807588, 42.62812185319396], [9.39000084802899, 43.009984849614824], [9.560016310269219, 42.152491970379458]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[143.910161981379559, 44.174099839853739], [144.613426548439662, 43.960882880217497], [145.320825230083159, 44.38473297787553], [145.543137241802782, 43.262088324550604], [144.059661899999895, 42.988358262700643], [143.183849725517263, 41.995214748699169], [141.611490920172486, 42.678790595056171], [141.067286411706647, 41.584593817707912], [139.955106235921079, 41.569555975910959], [139.817543573160009, 42.563758856774484], [140.312087030193283, 43.333272610032651], [141.380548944260084, 43.388824774746382], [141.67195234595394, 44.772125352551484], [141.96764489152801, 45.551483466161358], [143.142870314709825, 44.510358384776964], [143.910161981379559, 44.174099839853739]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-63.6645, 46.55001], [-62.9393, 46.41587], [-62.01208, 46.44314], [-62.50391, 46.033390000000111], [-62.87433, 45.96818], [-64.1428, 46.39265], [-64.39261, 46.72747], [-64.01486, 47.03601], [-63.6645, 46.55001]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-61.806305, 49.105060000000122], [-62.29318, 49.08717], [-63.58926, 49.40069], [-64.51912, 49.87304], [-64.17322, 49.95718], [-62.85829, 49.706410000000119], [-61.835585, 49.28855], [-61.806305, 49.105060000000122]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-123.510001587551173, 48.510010891303409], [-124.012890788399588, 48.370846259141501], [-125.655012777338371, 48.825004584338586], [-125.954994466792769, 49.179995835967645], [-126.850004435871895, 49.530000311880514], [-127.029993449544421, 49.814995835970052], [-128.059336304366241, 49.994959011426573], [-128.44458410710206, 50.539137681676038], [-128.358413656255436, 50.770648098343685], [-127.308581096029982, 50.55257355407204], [-126.695000977212317, 50.400903225295394], [-125.755006673823203, 50.295018215529467], [-125.415001587558805, 49.95000051533259], [-124.920768189119357, 49.47527497008349], [-123.922508708321118, 49.062483628935809], [-123.510001587551173, 48.510010891303409]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-56.134035814017039, 50.687009792679305], [-56.795881720595276, 49.81230866149096], [-56.143105027884218, 50.150117499382844], [-55.471492275602969, 49.935815334668462], [-55.822401089081012, 49.587128607778993], [-54.935142584845579, 49.31301097268684], [-54.473775397343701, 49.556691189159181], [-53.476549445191267, 49.249138902374057], [-53.786013759971325, 48.516780503933717], [-53.086133999226263, 48.687803656603535], [-52.958648240762244, 48.157164211614401], [-52.648098720904102, 47.53554840757559], [-53.069158291218343, 46.655498765645035], [-53.521456264852958, 46.61829173439483], [-54.178935512902541, 46.807065741556983], [-53.961868659060571, 47.625207017601923], [-54.240482143762051, 47.752279364607631], [-55.400773078011497, 46.884993801453135], [-55.997480841685842, 46.919720363953218], [-55.291219041552694, 47.389562486350997], [-56.250798712780636, 47.632545070987391], [-57.325229254777014, 47.572807115258087], [-59.26601518414688, 47.603347886742398], [-59.419494188053704, 47.899453843774864], [-58.796586473207498, 48.251525376979373], [-59.231624518456556, 48.523188381537892], [-58.391804979065228, 49.125580552764262], [-57.358689744685961, 50.718274034215938], [-56.738650071832012, 51.28743825947862], [-55.870976935435209, 51.632094224649279], [-55.406974249886616, 51.588272610065644], [-55.600218268442092, 51.317074693398013], [-56.134035814017039, 50.687009792679305]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-132.710007884431207, 54.040009315423447], [-132.710008504905716, 54.040009263721345], [-132.710007697461492, 54.040009263721345], [-132.710007884431207, 54.040009315423447], [-131.749989584003401, 54.120004380909137], [-132.049480347350993, 52.984621487024526], [-131.179042521826602, 52.180432847698285], [-131.577829549823122, 52.182370713909251], [-132.180428426778576, 52.639707139692405], [-132.549992432313871, 53.100014960332231], [-133.054611178755522, 53.411468817755292], [-133.239664482792705, 53.851080227262315], [-133.180004041711811, 54.169975490935315], [-132.710007884431207, 54.040009315423447]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[143.648007440362818, 50.747600409541604], [144.654147577085695, 48.976390692737482], [143.173927850517288, 49.306551418650372], [142.558668247650104, 47.861575018905], [143.533492466404113, 46.836728013692493], [143.505277134372619, 46.137907619809482], [142.747700636973917, 46.740764878926456], [142.092030064054569, 45.966755276058791], [141.906925083585037, 46.805928860046549], [142.018442824470952, 47.780132961612935], [141.904444614835057, 48.85918854429957], [142.135800002205741, 49.615163072297349], [142.179983351815309, 50.952342434281917], [141.594075962490052, 51.935434882202543], [141.68254601457366, 53.301966457728867], [142.606934035410831, 53.762145087287905], [142.209748976815462, 54.225475979216952], [142.65478641171299, 54.365880845753878], [142.91461551327663, 53.704577541714826], [143.260847609632123, 52.740760403039133], [143.235267775647657, 51.756660264688747], [143.648007440362818, 50.747600409541604]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-6.788856573910806, 52.260117906292436], [-8.561616583683502, 51.669301255899398], [-9.97708574059024, 51.820454820353035], [-9.166282517930767, 52.864628811242653], [-9.688524542672383, 53.881362616585363], [-8.327987433291923, 54.664518947968588], [-7.572167934590993, 55.131622219454982], [-6.73384701173606, 55.172860012423797], [-5.661948614921926, 54.554603176483766], [-6.197884894220977, 53.867565009163343], [-6.032985398777498, 53.153164170944308], [-6.788856573910806, 52.260117906292436]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[12.690006137755631, 55.609990953180869], [12.089991082414713, 54.800014553438018], [11.043543328504313, 55.364863796604169], [10.903913608451688, 55.779954738988835], [12.370904168353377, 56.111407375708922], [12.690006137755631, 55.609990953180869]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-153.006314053336922, 57.115842190165893], [-154.005090298458128, 56.734676825581062], [-154.516402757770095, 56.992748928446701], [-154.670992804971149, 57.4611957871725], [-153.762779507441508, 57.816574612043752], [-153.228729417921187, 57.968968410872435], [-152.564790615835136, 57.901427313866975], [-152.141147223906444, 57.591058661522084], [-153.006314053336922, 57.115842190165893]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-3.005004848635195, 58.635000108466244], [-4.073828497728016, 57.553024807355229], [-3.055001796877661, 57.69001902936094], [-1.959280564776833, 57.68479970969949], [-2.219988165689387, 56.870017401753529], [-3.119003058271119, 55.973793036515474], [-2.085009324543051, 55.909998480851186], [-1.114991013992238, 54.624986477265395], [-0.430484991854115, 54.464376125702159], [0.184981316742125, 53.325014146530947], [0.469976840831777, 52.929999498091888], [1.681530795914824, 52.739520168664001], [1.559987827164292, 52.099998480836092], [1.050561557630999, 51.806760565795685], [1.449865349950301, 51.289427802121878], [0.550333693045673, 50.765738837276075], [-0.78751746255864, 50.774988918656305], [-2.489997524414463, 50.500018622431242], [-2.956273972983951, 50.696879991247101], [-3.617448085942328, 50.228355617872722], [-4.542507900399158, 50.341837063185665], [-5.245023159191049, 49.959999904981174], [-5.776566941745301, 50.159677639356914], [-4.309989793301924, 51.210001125689161], [-3.414850633142038, 51.426008612669165], [-4.984367234710788, 51.593466091510976], [-5.267295701508885, 51.99140045837467], [-4.222346564134853, 52.301355699261364], [-4.770013393564113, 52.84000499125554], [-4.579999152026915, 53.495003770555087], [-3.092079637047021, 53.404440822963551], [-2.945008510744373, 53.984999701546599], [-3.630005458989331, 54.615012925833014], [-4.844169073903004, 54.790971177786929], [-5.082526617849226, 55.061600653699372], [-4.719112107756644, 55.508472601943282], [-5.047980922862138, 55.783985500707445], [-5.58639767091114, 55.311146145236819], [-5.644998745130181, 56.275014960344805], [-6.149980841486354, 56.785009670633457], [-5.786824713555205, 57.818848375064562], [-5.009998745127575, 58.630013332750053], [-4.211494513353585, 58.550845038479054], [-3.005004848635195, 58.635000108466244]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-165.579164191733611, 59.90998688418756], [-166.192770148767295, 59.754440822989068], [-166.848337368822001, 59.941406155021042], [-167.455277066090076, 60.213069159579476], [-166.467792121424651, 60.384169826897875], [-165.674429694663672, 60.293606879306338], [-165.579164191733611, 59.90998688418756]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-79.26582, 62.158675], [-79.65752, 61.63308], [-80.09956, 61.718100000000135], [-80.36215, 62.01649], [-80.315395, 62.085565000000116], [-79.92939, 62.3856], [-79.52002, 62.36371], [-79.26582, 62.158675]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-81.898249999999877, 62.7108], [-83.06857, 62.159220000000118], [-83.774619999999885, 62.18231], [-83.99367, 62.4528], [-83.25048, 62.91409000000013], [-81.87699, 62.90458], [-81.898249999999877, 62.7108]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-171.731656867539414, 63.78251536727592], [-171.114433560245345, 63.59219106714491], [-170.491112433940742, 63.694975490973604], [-169.682505459653584, 63.431115627691241], [-168.68943946030069, 63.297506212000513], [-168.771940884454722, 63.188598130945451], [-169.529439867205156, 62.976931464277897], [-170.290556200215974, 63.194437567794552], [-170.67138566799099, 63.375821845138972], [-171.553063117538784, 63.317789211675091], [-171.79111060289128, 63.405845852300587], [-171.731656867539414, 63.78251536727592]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-85.161307949549951, 65.657284654392896], [-84.975763719405876, 65.21751821558891], [-84.464012010419424, 65.371772365980263], [-83.882626308919669, 65.109617824963635], [-82.787576870438699, 64.766693020274687], [-81.642013719392651, 64.455135809986956], [-81.553440314444373, 63.979609280037238], [-80.817361212878865, 64.057485663501097], [-80.103451300766523, 63.725981350348604], [-80.991019863595596, 63.411246039474889], [-82.547178107417011, 63.651722317145243], [-83.108797573565056, 64.101875718839807], [-84.100416632813875, 63.569711819098103], [-85.523404710619047, 63.052379055424097], [-85.866768764982396, 63.637252916103563], [-87.22198320183665, 63.541238104905148], [-86.352759772471273, 64.035833238370714], [-86.224886440765147, 64.822916978608276], [-85.883847825854787, 65.738778388117055], [-85.161307949549951, 65.657284654392896]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-14.508695441129134, 66.455892239031385], [-14.739637417041536, 65.808748277440372], [-13.609732224979695, 65.126671047619936], [-14.909833746794817, 64.364081936288642], [-17.794438035543351, 63.678749091233925], [-18.656245896874964, 63.496382961675835], [-19.972754685942732, 63.643634955491535], [-22.76297197111009, 63.960178941495371], [-21.778484259517597, 64.402115790455468], [-23.955043911219065, 64.891129869233481], [-22.184402635170301, 65.084968166760319], [-22.227423265053233, 65.378593655042749], [-24.326184047939222, 65.611189276788451], [-23.650514695723075, 66.262519029395236], [-22.134922451250816, 66.410468655046799], [-20.576283738679479, 65.732112128351531], [-19.056841600001576, 66.276600857194893], [-17.798623826559009, 65.993853257909763], [-16.16781897629204, 66.526792304135796], [-14.508695441129134, 66.455892239031385]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-75.86588, 67.14886], [-76.98687, 67.09873], [-77.2364, 67.588090000000108], [-76.81166, 68.148560000000117], [-75.89521, 68.28721], [-75.1145, 68.01036], [-75.10333, 67.58202], [-75.21597, 67.44425], [-75.86588, 67.14886]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-175.01425, 66.58435], [-174.33983, 66.3355600000001], [-174.57182, 67.06219], [-171.85731, 66.91308], [-169.89958, 65.97724], [-170.89107, 65.54139], [-172.53025, 65.43791], [-172.555, 64.46079], [-172.95533, 64.25269], [-173.89184, 64.2826], [-174.65392, 64.63125], [-175.98353, 64.92288], [-176.20716, 65.35667], [-177.22266, 65.52024], [-178.35993, 65.39052], [-178.90332, 65.74044], [-178.68611, 66.11211], [-179.88377, 65.87456], [-179.43268, 65.40411], [-180.0, 64.979708702198451], [-180.0, 68.963636363636454], [-177.55, 68.200000000000131], [-174.92825, 67.20589], [-175.01425, 66.58435]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-95.647681203800431, 69.107690358321776], [-96.269521203800593, 68.757040358321746], [-97.617401203800569, 69.060030358321711], [-98.431801203800546, 68.950700358321711], [-99.797401203800547, 69.400030358321715], [-98.917401203800551, 69.710030358321887], [-98.218261203800608, 70.143540358321843], [-97.157401203800475, 69.860030358321893], [-96.557401203800453, 69.680030358321773], [-96.257401203800441, 69.490030358321775], [-95.647681203800431, 69.107690358321776]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[180.000000000000142, 70.832199208546683], [178.903425000000112, 70.78114], [178.7253, 71.0988], [180.000000000000142, 71.51571433642826], [180.000000000000142, 70.832199208546683]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-178.69378, 70.89302], [-180.0, 70.832199208546683], [-180.0, 71.51571433642826], [-179.871875, 71.55762], [-179.02433, 71.55553], [-177.577945, 71.26948], [-177.663575, 71.13277], [-178.69378, 70.89302]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-90.547119140624943, 69.497680664062585], [-90.551513671874943, 68.475097656250028], [-89.215087890624915, 69.258728027343864], [-88.019592285156193, 68.615112304687557], [-88.317504882812443, 67.87347412109375], [-87.350097656249943, 67.198730468750057], [-86.306091308593665, 67.921508789062642], [-85.576599121093722, 68.784484863281307], [-85.521911621093778, 69.882080078124972], [-84.100769042968722, 69.805480957031278], [-82.622497558593722, 69.658325195312614], [-81.280395507812472, 69.162109375000028], [-81.220214843749972, 68.665710449218835], [-81.964294433593665, 68.13250732421875], [-81.259277343750028, 67.5972900390625], [-81.386474609374915, 67.110900878906222], [-83.344482421875028, 66.41168212890625], [-84.735412597656193, 66.257324218750114], [-85.769409179687443, 66.558288574218864], [-86.067626953124972, 66.056274414062614], [-87.031372070312415, 65.213073730468892], [-87.323181152343665, 64.775695800781364], [-88.482910156249972, 64.099121093750085], [-89.914428710937415, 64.032714843750057], [-90.703979492187443, 63.610290527343864], [-90.770019531249972, 62.960327148437472], [-91.933410644531222, 62.835083007812557], [-93.156982421874972, 62.024719238281364], [-94.241516113281165, 60.898681640625], [-94.629272460937415, 60.110290527343835], [-94.684509277343665, 58.948913574218835], [-93.215026855468693, 58.7821044921875], [-92.764587402343693, 57.845703125000085], [-92.296997070312443, 57.087097167968807], [-90.897705078124972, 57.284729003906335], [-89.039489746093722, 56.8516845703125], [-88.039794921874943, 56.471679687500142], [-87.324218749999886, 55.99908447265625], [-86.071228027343778, 55.723876953125028], [-85.011779785156165, 55.302673339843892], [-83.360473632812415, 55.244873046875142], [-82.272827148437415, 55.148315429687585], [-82.436218261718693, 54.28228759765625], [-82.125, 53.277099609375], [-81.400695800781193, 52.157897949218778], [-79.912902832031222, 51.208496093750085], [-79.143005371093722, 51.533874511718835], [-78.601928710937415, 52.562072753906278], [-79.124206542968693, 54.141479492187642], [-79.829589843749915, 54.667724609375142], [-78.228698730468693, 55.136474609375085], [-77.095581054687472, 55.8375244140625], [-76.541381835937472, 56.5343017578125], [-76.623107910156222, 57.202697753906307], [-77.302185058593693, 58.052124023437642], [-78.516906738281193, 58.804687500000057], [-77.336669921874972, 59.852722167968778], [-77.772705078124915, 60.757873535156364], [-78.106811523437528, 62.319702148437528], [-77.41058349609375, 62.550476074218778], [-75.696228027343722, 62.278503417968778], [-74.668212890624915, 62.181091308593835], [-73.839904785156222, 62.443908691406307], [-72.908508300781165, 62.1051025390625], [-71.677001953124915, 61.525329589843807], [-71.373718261718722, 61.137084960937472], [-69.590393066406193, 61.061523437500114], [-69.620300292968722, 60.221313476562614], [-69.287902832031193, 58.957275390625057], [-68.374511718749972, 58.801086425781335], [-67.649780273437443, 58.212097167968835], [-66.201782226562415, 58.767272949218778], [-65.245178222656278, 59.870727539062585], [-64.583496093749915, 60.335693359375057], [-63.804687499999972, 59.442687988281307], [-62.502380371093665, 58.167114257812614], [-61.396484375, 56.967529296875085], [-61.798583984374915, 56.339477539062614], [-60.468505859374943, 55.775512695312642], [-59.569580078124943, 55.2041015625], [-57.975097656249915, 54.945495605468778], [-57.333190917968665, 54.626525878906392], [-56.936889648437472, 53.7802734375], [-56.158081054687443, 53.64752197265625], [-55.756286621093665, 53.270507812500057], [-55.683288574218722, 52.146728515625085], [-56.409179687499915, 51.77069091796875], [-57.126892089843665, 51.419677734375114], [-58.774780273437443, 51.064270019531307], [-60.033081054687472, 50.242919921875028], [-61.723571777343722, 50.080505371093864], [-63.862487792968722, 50.291076660156364], [-65.363281249999972, 50.298278808593807], [-66.398986816406165, 50.228881835937585], [-67.236328124999915, 49.511474609375028], [-68.511108398437415, 49.068481445312614], [-69.95361328125, 47.744873046875028], [-71.104492187499972, 46.821716308593778], [-70.255187988281165, 46.986083984374972], [-68.650024414062472, 48.300109863281278], [-66.55242919921875, 49.133117675781335], [-65.056213378906165, 49.232910156250142], [-64.171020507812415, 48.742492675781335], [-65.115478515625028, 48.070922851562614], [-64.798522949218722, 46.993103027343864], [-64.472106933593665, 46.238525390625085], [-63.173278808593665, 45.73907470703125], [-61.520690917968693, 45.883911132812557], [-60.518127441406222, 47.007873535156278], [-60.448608398437472, 46.282714843750028], [-59.802795410156222, 45.92047119140625], [-61.039794921874972, 45.265319824218778], [-63.254699707031222, 44.670288085937585], [-64.246582031249943, 44.265502929687528], [-65.364074707031165, 43.545288085937528], [-66.123413085937443, 43.618713378906278], [-66.161682128906222, 44.465087890625114], [-64.425476074218665, 45.292114257812614], [-66.0260009765625, 45.259277343750085], [-67.137390136718722, 45.13751220703125], [-66.964599609374972, 44.809692382812557], [-68.032470703124943, 44.325317382812642], [-69.059997558593722, 43.980102539062472], [-70.116088867187472, 43.684082031249972], [-70.690002441406193, 43.030090332031278], [-70.814880371093665, 42.865295410156278], [-70.825012207031165, 42.335083007812557], [-70.494995117187472, 41.805114746093864], [-70.080017089843693, 41.780090332031335], [-70.184997558593665, 42.145080566406278], [-69.884887695312443, 41.922912597656307], [-69.965026855468722, 41.637084960937642], [-70.640014648437443, 41.475097656250085], [-71.120300292968693, 41.494506835937528], [-71.859985351562472, 41.32012939453125], [-72.294982910156165, 41.270080566406364], [-72.876403808593722, 41.220703125000114], [-73.710021972656222, 40.931091308593807], [-72.241210937499915, 41.119506835937528], [-71.945007324218665, 40.930114746093778], [-73.344909667968665, 40.630126953125057], [-73.981994628906193, 40.628112792968892], [-73.952270507812472, 40.750671386718807], [-74.256713867187415, 40.473510742187557], [-73.962402343749915, 40.42767333984375], [-74.178405761718665, 39.709289550781335], [-74.906005859374943, 38.939514160156364], [-74.980407714843722, 39.196472167968722], [-75.200012207031278, 39.248474121093835], [-75.528076171875028, 39.498474121093807], [-75.320007324218722, 38.960083007812585], [-75.083496093749972, 38.781311035156307], [-75.056701660156222, 38.404113769531364], [-75.377380371093665, 38.015502929687528], [-75.940185546874972, 37.216918945312557], [-76.03118896484375, 37.2567138671875], [-75.721984863281222, 37.93707275390625], [-76.232788085937472, 38.319274902343892], [-76.349975585937472, 39.150085449218864], [-76.542724609374972, 38.717712402343892], [-76.329284667968722, 38.083312988281392], [-76.960021972656193, 38.23291015625], [-76.301574707031165, 37.918090820312585], [-76.258728027343722, 36.966491699218807], [-75.971801757812415, 36.897277832031335], [-75.867980957031193, 36.551330566406392], [-75.727478027343693, 35.550720214843864], [-76.363098144531165, 34.808471679687528], [-77.397583007812472, 34.512084960937528], [-78.054870605468722, 33.925476074218864], [-78.554321289062472, 33.861328125000057], [-79.060607910156222, 33.49407958984375], [-79.2034912109375, 33.158508300781364], [-80.301269531249915, 32.509277343750142], [-80.864990234374915, 32.033325195312557], [-81.336303710937415, 31.440490722656278], [-81.490417480468665, 30.730102539062642], [-81.313720703124943, 30.035522460937528], [-80.979980468749972, 29.180114746093807], [-80.535583496093636, 28.472106933593807], [-80.530029296874972, 28.040100097656307], [-80.056518554687443, 26.880126953125142], [-80.088012695312472, 26.205688476562642], [-80.131469726562443, 25.816894531250057], [-80.380981445312472, 25.206298828125028], [-80.679992675781165, 25.080078125000057], [-81.172119140624943, 25.201293945312642], [-81.330017089843665, 25.640075683593807], [-81.710021972656222, 25.870117187500085], [-82.239990234374886, 26.730102539062528], [-82.705078124999972, 27.495117187500028], [-82.855285644531222, 27.88629150390625], [-82.650024414062443, 28.550109863281307], [-82.929992675781165, 29.100097656250142], [-83.709594726562472, 29.93670654296875], [-84.099975585937443, 30.090087890625028], [-85.108825683593722, 29.63629150390625], [-85.287780761718693, 29.686096191406278], [-85.773010253906193, 30.152709960937585], [-86.400024414062415, 30.400085449218835], [-87.530273437499943, 30.27447509765625], [-88.417785644531222, 30.384887695312642], [-89.180480957031165, 30.316101074218778], [-89.604919433593722, 30.176330566406307], [-89.413696289062472, 29.894287109375114], [-89.429992675781193, 29.488708496093807], [-89.217590332031278, 29.291076660156278], [-89.408203124999915, 29.159729003906335], [-89.779296874999943, 29.307128906250057], [-90.154602050781222, 29.11749267578125], [-90.880187988281165, 29.148681640625085], [-91.626708984374972, 29.677124023437528], [-92.499084472656193, 29.552307128906364], [-93.226379394531193, 29.783874511718892], [-93.848388671874972, 29.713684082031392], [-94.690002441406222, 29.480102539062528], [-95.600280761718722, 28.738708496093778], [-96.593994140624943, 28.307495117187557], [-97.140014648437415, 27.830078125000057], [-97.369995117187443, 27.380126953125028], [-97.380004882812443, 26.690124511718864], [-97.330017089843693, 26.210083007812642], [-97.140197753906136, 25.869506835937528], [-97.138610839843693, 25.867919921874972], [-97.141784667968693, 25.865905761718807], [-97.528076171874915, 24.992126464843778], [-97.702880859374915, 24.272277832031278], [-97.776000976562472, 22.932678222656392], [-97.872375488281278, 22.444274902343778], [-97.698974609374915, 21.898681640625114], [-97.388977050781193, 21.411071777343778], [-97.189270019531193, 20.635498046875], [-96.525512695312443, 19.890930175781307], [-96.292114257812415, 19.320495605468835], [-95.900878906249972, 18.828125], [-94.838989257812415, 18.562683105468835], [-94.425720214843693, 18.144287109375142], [-93.548583984374915, 18.423889160156392], [-92.786071777343722, 18.524902343750057], [-92.037292480468665, 18.7047119140625], [-91.407897949218778, 18.876098632812614], [-90.771789550781193, 19.284118652343892], [-90.533508300781222, 19.867492675781307], [-90.451477050781222, 20.707519531250085], [-90.278625488281165, 20.999877929687557], [-89.601318359374915, 21.261718750000028], [-88.543884277343665, 21.493713378906335], [-87.658386230468636, 21.45892333984375], [-87.051879882812415, 21.543518066406364], [-86.812011718749915, 21.331481933593722], [-86.845886230468665, 20.849914550781278], [-87.383300781249943, 20.255493164062614], [-87.620971679687472, 19.646484375000114], [-87.436706542968693, 19.472473144531307], [-87.586486816406165, 19.040100097656335], [-87.837219238281165, 18.259887695312472], [-88.090576171874972, 18.516723632812585], [-88.299987792968722, 18.500122070312585], [-88.296325683593665, 18.353271484375085], [-88.106811523437415, 18.348693847656335], [-88.123413085937415, 18.076721191406364], [-88.285278320312528, 17.644287109375057], [-88.197875976562472, 17.489501953125142], [-88.302612304687415, 17.1317138671875], [-88.239501953125028, 17.036071777343864], [-88.355407714843722, 16.5308837890625], [-88.551818847656165, 16.265502929687614], [-88.732421874999915, 16.233703613281392], [-88.930603027343665, 15.887329101562642], [-88.604614257812472, 15.706481933593864], [-88.518310546874972, 15.855529785156278], [-88.224975585937415, 15.72772216796875], [-88.121093749999972, 15.688720703125114], [-87.901794433593665, 15.864501953125028], [-87.615600585937443, 15.878906250000114], [-87.522888183593665, 15.797302246093892], [-87.367675781249972, 15.846923828125], [-86.903198242187415, 15.7567138671875], [-86.440917968749972, 15.782897949218807], [-86.119201660156278, 15.893493652343778], [-86.001892089843693, 16.005493164062557], [-85.683288574218693, 15.953674316406364], [-85.443969726562415, 15.885681152343778], [-85.182373046874915, 15.9093017578125], [-84.983703613281222, 15.995910644531278], [-84.526977539062472, 15.857299804687585], [-84.368225097656165, 15.8350830078125], [-84.062988281249943, 15.648315429687642], [-83.773986816406165, 15.424072265625], [-83.410400390624943, 15.270874023437642], [-83.147216796874972, 14.995910644531307], [-83.233215332031222, 14.899902343749972], [-83.284179687499915, 14.676696777343864], [-83.182128906249915, 14.310729980468736], [-83.412475585937443, 13.9700927734375], [-83.519775390624943, 13.567687988281293], [-83.552185058593693, 13.127075195312528], [-83.498474121093665, 12.869323730468892], [-83.473327636718778, 12.419128417968835], [-83.626098632812443, 12.320922851562614], [-83.719604492187472, 11.893127441406378], [-83.650878906249915, 11.629089355468821], [-83.855407714843722, 11.373291015625057], [-83.808898925781165, 11.103088378906364], [-83.655578613281193, 10.938903808593906], [-83.402282714843693, 10.395507812499986], [-83.015686035156222, 9.993103027343778], [-82.546203613281222, 9.566284179687571], [-82.187072753906165, 9.207519531250071], [-82.207580566406165, 8.9957275390625], [-81.808593750000028, 8.950683593750156], [-81.714111328124972, 9.032104492187628], [-81.439208984374915, 8.786315917968764], [-80.947326660156278, 8.858520507812614], [-80.521911621093693, 9.111083984375142], [-79.914611816406193, 9.312683105468892], [-79.573303222656222, 9.611694335937585], [-79.021179199218665, 9.552917480468807], [-79.058410644531222, 9.454711914062571], [-78.500915527343722, 9.420471191406236], [-78.055908203124972, 9.2476806640625], [-77.7294921875, 8.9468994140625], [-77.353271484374915, 8.670471191406378], [-76.836608886718778, 8.638671875], [-76.086303710937472, 9.336914062500142], [-75.674621582031165, 9.443298339843849], [-75.664672851562472, 9.77410888671875], [-75.480407714843693, 10.619079589843764], [-74.906921386718665, 11.0831298828125], [-74.276672363281193, 11.102111816406349], [-74.197204589843665, 11.31048583984375], [-73.414672851562472, 11.227111816406236], [-72.627807617187472, 11.732116699218849], [-72.238220214843722, 11.955688476562628], [-71.754089355468693, 12.437316894531406], [-71.399780273437415, 12.376098632812486], [-71.137390136718722, 12.11309814453125], [-71.331604003906165, 11.776306152343793], [-71.359985351562472, 11.540100097656406], [-71.947021484374943, 11.423278808593821], [-71.620788574218693, 10.969482421875028], [-71.632995605468693, 10.446472167968764], [-72.074096679687472, 9.865722656250057], [-71.695617675781165, 9.072326660156364], [-71.264587402343693, 9.137329101562571], [-71.039978027343693, 9.860107421875156], [-71.350097656249915, 10.211914062500057], [-71.400573730468778, 10.969116210937543], [-70.155212402343722, 11.375488281250142], [-70.293823242187443, 11.846923828125085], [-69.943176269531193, 12.162292480468736], [-69.584289550781193, 11.459716796875071], [-68.882995605468722, 11.443481445312557], [-68.233276367187472, 10.885681152343892], [-68.194091796874915, 10.554687500000085], [-67.296203613281165, 10.545898437500071], [-66.227783203124972, 10.648681640625043], [-65.65521240234375, 10.200927734375128], [-64.890380859374915, 10.077270507812571], [-64.329406738281222, 10.389709472656335], [-64.317993164062415, 10.641479492187585], [-63.079284667968693, 10.701721191406307], [-61.880920410156222, 10.715698242187614], [-62.7301025390625, 10.420288085937642], [-62.388488769531222, 9.948303222656293], [-61.588684082031193, 9.87310791015625], [-60.830505371093665, 9.381286621093835], [-60.671203613281222, 8.580322265625014], [-60.150085449218693, 8.602905273437585], [-59.758300781249943, 8.367126464843807], [-59.101684570312415, 7.99932861328125], [-58.482910156249972, 7.347717285156378], [-58.454895019531165, 6.832885742187585], [-58.078125, 6.809082031250128], [-57.542175292968693, 6.321289062500043], [-57.147399902343665, 5.973083496093864], [-55.949279785156222, 5.772888183593764], [-55.841796874999972, 5.953125], [-55.033203124999915, 6.025329589843849], [-53.958007812499972, 5.756530761718807], [-53.618408203125, 5.646484375000057], [-52.882080078125, 5.409912109375], [-51.823303222656222, 4.565917968750014], [-51.657775878906165, 4.156311035156349], [-51.317077636718722, 4.203491210937656], [-51.069702148437472, 3.650512695312628], [-50.508789062499915, 1.901489257812642], [-49.973999023437443, 1.736511230468892], [-49.947082519531193, 1.046325683593793], [-50.699279785156222, 0.223083496093764], [-50.388183593749943, -0.078369140624986], [-48.620483398437472, -0.235412597656207], [-48.584411621093722, -1.237792968749929], [-47.824890136718693, -0.581604003906207], [-46.566589355468665, -0.940979003906236], [-44.905700683593693, -1.551696777343722], [-44.417602539062443, -2.137695312499943], [-44.581604003906165, -2.691284179687429], [-43.418701171874915, -2.383117675781236], [-41.472595214843722, -2.911987304687443], [-39.97857666015625, -2.872985839843707], [-38.500305175781165, -3.700622558593736], [-37.223205566406222, -4.820922851562415], [-36.452880859374943, -5.109375], [-35.5977783203125, -5.149475097656179], [-35.235412597656165, -5.464904785156222], [-34.895996093749915, -6.738220214843693], [-34.729980468749943, -7.343200683593665], [-35.128173828124972, -8.996398925781165], [-35.636901855468665, -9.649291992187429], [-37.0465087890625, -11.040710449218679], [-37.683593749999943, -12.171203613281193], [-38.423889160156222, -13.038085937499929], [-38.673889160156193, -13.057678222656222], [-38.95318603515625, -13.793395996093679], [-38.882324218749972, -15.666992187499943], [-39.161010742187415, -17.208374023437486], [-39.267272949218693, -17.867675781249929], [-39.583496093749915, -18.262207031249929], [-39.760803222656193, -19.599121093749972], [-40.774719238281165, -20.904479980468778], [-40.944702148437443, -21.937316894531207], [-41.754089355468693, -22.370605468749986], [-41.988281249999915, -22.970092773437486], [-43.074707031249915, -22.967712402343651], [-44.647827148437443, -23.351989746093778], [-45.352111816406193, -23.796813964843722], [-46.472106933593693, -24.088989257812528], [-47.64898681640625, -24.885192871093679], [-48.495483398437415, -25.877014160156179], [-48.6409912109375, -26.623718261718651], [-48.47467041015625, -27.175903320312415], [-48.6614990234375, -28.186096191406222], [-48.888427734374943, -28.674072265624957], [-49.587280273437443, -29.224487304687415], [-50.696899414062472, -30.984375], [-51.576171874999943, -31.777709960937457], [-52.256103515624915, -32.245300292968707], [-52.712097167968693, -33.196594238281207], [-53.373596191406222, -33.768310546875043], [-53.806396484374915, -34.396789550781207], [-54.935791015624972, -34.952575683593778], [-55.674011230468693, -34.752685546875028], [-56.215209960937443, -34.859802246093707], [-57.139709472656278, -34.430480957031236], [-57.817810058593722, -34.462524414062486], [-58.427001953124972, -33.909484863281222], [-58.495422363281165, -34.431518554687472], [-57.225769042968722, -35.288024902343736], [-57.362304687499943, -35.977416992187457], [-56.737487792968722, -36.413085937499901], [-56.788208007812472, -36.901489257812429], [-57.749084472656165, -38.183898925781179], [-59.231811523437472, -38.720214843749901], [-61.237426757812415, -38.928405761718651], [-62.335876464843665, -38.827697753906151], [-62.125793457031193, -39.424072265625028], [-62.330505371093722, -40.172607421875028], [-62.145996093749915, -40.676879882812415], [-62.745788574218778, -41.028686523437415], [-63.770507812499943, -41.166809082031193], [-64.732116699218693, -40.802612304687543], [-65.117980957031222, -41.064270019531165], [-64.978576660156222, -42.057983398437401], [-64.303405761718665, -42.359008789062543], [-63.755920410156278, -42.043701171874929], [-63.458007812499943, -42.563110351562457], [-64.378784179687472, -42.873474121093665], [-65.181823730468665, -43.495300292968693], [-65.328796386718665, -44.501281738281236], [-65.565185546875, -45.036804199218679], [-66.509887695312472, -45.039611816406207], [-67.293823242187415, -45.551879882812415], [-67.580505371093665, -46.301696777343707], [-66.596984863281165, -47.033874511718651], [-65.640991210937472, -47.236083984374929], [-65.985107421874915, -48.133300781249929], [-67.166198730468722, -48.697326660156207], [-67.816101074218693, -49.869689941406214], [-68.728698730468722, -50.264221191406214], [-69.138488769531222, -50.732482910156222], [-68.815490722656222, -51.771118164062457], [-68.1500244140625, -52.349975585937528], [-68.571472167968665, -52.299377441406229], [-69.461303710937472, -52.291870117187415], [-69.942687988281193, -52.537902832031236], [-70.845092773437472, -52.899169921874915], [-71.006286621093665, -53.833190917968658], [-71.429809570312415, -53.856384277343672], [-72.557922363281165, -53.531372070312422], [-73.70269775390625, -52.835083007812415], [-74.946777343749915, -52.262695312499922], [-75.260009765624943, -51.629272460937429], [-74.976623535156165, -51.043395996093736], [-75.479675292968722, -50.378295898437479], [-75.607971191406165, -48.673706054687486], [-75.182678222656222, -47.711914062499972], [-74.126586914062472, -46.939208984374943], [-75.644409179687415, -46.647583007812528], [-74.692077636718693, -45.763977050781165], [-74.351684570312415, -44.103027343749986], [-73.240295410156193, -44.454895019531278], [-72.717712402343722, -42.383300781249929], [-73.388916015624972, -42.117492675781151], [-73.701293945312415, -43.365783691406179], [-74.331909179687472, -43.224975585937472], [-74.017883300781222, -41.794799804687429], [-73.677124023437472, -39.942199707031278], [-73.217590332031165, -39.258605957031278], [-73.505493164062415, -38.282897949218693], [-73.588012695312443, -37.156311035156179], [-73.166687011718693, -37.123779296874915], [-72.553100585937443, -35.508789062499972], [-71.861694335937472, -33.909118652343651], [-71.438476562499972, -32.418884277343722], [-71.668701171875, -30.920593261718651], [-71.369995117187472, -30.095703125000028], [-71.489807128906165, -28.861389160156222], [-70.905090332031222, -27.640380859375043], [-70.724975585937415, -25.705871582031222], [-70.403991699218693, -23.628906249999915], [-70.091186523437443, -21.393310546874929], [-70.164428710937443, -19.756408691406165], [-70.372497558593665, -18.347900390624972], [-71.375183105468665, -17.773803710937486], [-71.461975097656165, -17.363403320312457], [-73.444519042968722, -16.359375], [-75.237792968749943, -15.265686035156193], [-76.009216308593693, -14.649291992187415], [-76.423400878906193, -13.823181152343707], [-76.259216308593665, -13.534973144531207], [-77.106201171874915, -12.222717285156207], [-78.092102050781165, -10.377685546874901], [-79.036926269531278, -8.386596679687486], [-79.445922851562415, -7.930786132812443], [-79.760498046874943, -7.194274902343707], [-80.537475585937472, -6.541687011718707], [-81.25, -6.136779785156236], [-80.926269531249943, -5.690490722656179], [-81.410888671874915, -4.736694335937486], [-81.099609374999915, -4.036376953124957], [-80.302490234374972, -3.404785156249986], [-79.770202636718693, -2.657470703124972], [-79.986511230468693, -2.220703124999986], [-80.368713378906165, -2.685180664062429], [-80.967712402343665, -2.246887207031193], [-80.764770507812415, -1.965026855468736], [-80.933593750000028, -1.057373046874943], [-80.583312988281278, -0.906677246093679], [-80.399291992187472, -0.283691406249972], [-80.020812988281165, 0.360473632812571], [-80.090576171874972, 0.768493652343864], [-79.542785644531222, 0.982910156250142], [-78.855285644531222, 1.380920410156335], [-78.990905761718722, 1.691284179687628], [-78.617797851562443, 1.766479492187486], [-78.662109374999915, 2.267272949218849], [-78.427612304687415, 2.629699707031378], [-77.931518554687443, 2.696716308593764], [-77.510375976562443, 3.325073242187585], [-77.127685546874943, 3.849670410156406], [-77.496276855468693, 4.087707519531406], [-77.307617187499915, 4.668090820312628], [-77.533203124999972, 5.582885742187656], [-77.318786621093778, 5.845275878906349], [-77.476684570312415, 6.691101074218778], [-77.881591796874972, 7.223876953125071], [-78.214904785156222, 7.512329101562628], [-78.429077148437472, 8.052124023437642], [-78.182006835937472, 8.319274902343835], [-78.435485839843722, 8.387695312500028], [-78.622070312499943, 8.718078613281307], [-79.120300292968722, 8.996093749999986], [-79.557800292968665, 8.932495117187585], [-79.760498046874943, 8.584472656250156], [-80.164489746093722, 8.333312988281264], [-80.382629394531165, 8.298522949218864], [-80.480712890624943, 8.090270996093892], [-80.003601074218722, 7.547485351562514], [-80.276611328125028, 7.419677734375], [-80.421081542968693, 7.271484375], [-80.886413574218722, 7.220520019531392], [-81.059509277343722, 7.817871093750114], [-81.189697265624972, 7.647888183593835], [-81.519470214843722, 7.706726074218736], [-81.721313476562415, 8.108886718750071], [-82.131408691406278, 8.175476074218849], [-82.390869140624915, 8.29248046875], [-82.820007324218722, 8.290893554687614], [-82.850891113281222, 8.073913574218764], [-82.965698242187443, 8.225097656250128], [-83.508422851562528, 8.446899414062599], [-83.711486816406222, 8.656921386718878], [-83.596313476562415, 8.830505371093892], [-83.632629394531222, 9.051513671875085], [-83.909912109374943, 9.290893554687599], [-84.303405761718778, 9.487487792968793], [-84.647583007812528, 9.615478515625057], [-84.713378906249943, 9.908081054687571], [-84.975585937499972, 10.08673095703125], [-84.911376953124943, 9.796081542968778], [-85.110900878906222, 9.557128906250071], [-85.339477539062528, 9.834472656250085], [-85.660705566406222, 9.933288574218849], [-85.797424316406193, 10.134887695312585], [-85.791687011718665, 10.439270019531307], [-85.659301757812415, 10.754272460937656], [-85.941711425781193, 10.895324707031321], [-85.712524414062415, 11.088500976562528], [-86.058410644531165, 11.403503417968878], [-86.525878906249972, 11.806884765625099], [-86.745910644531136, 12.144104003906349], [-87.167480468750028, 12.458312988281236], [-87.668518066406165, 12.909912109375114], [-87.557495117187415, 13.064697265625028], [-87.392395019531222, 12.914123535156378], [-87.316589355468636, 12.984680175781364], [-87.489379882812472, 13.297485351562614], [-87.793090820312528, 13.384521484374986], [-87.904113769531193, 13.149108886718878], [-88.483276367187443, 13.163879394531278], [-88.843200683593722, 13.259704589843849], [-89.256713867187472, 13.458679199218878], [-89.812377929687415, 13.520690917968906], [-90.095581054687472, 13.735473632812486], [-90.608581542968665, 13.909912109375085], [-91.232421875000028, 13.927917480468736], [-91.689697265624943, 14.126281738281236], [-92.227722167968722, 14.538879394531264], [-93.359375, 15.615478515625114], [-93.875183105468722, 15.940307617187614], [-94.691589355468778, 16.20111083984375], [-95.250183105468722, 16.128295898437528], [-96.053405761718665, 15.752075195312642], [-96.557373046874943, 15.65350341796875], [-97.263610839843665, 15.9171142578125], [-98.013000488281136, 16.107299804687528], [-98.947692871093736, 16.566101074218864], [-99.697387695312401, 16.706298828125114], [-100.829528808593778, 17.171081542968835], [-101.666076660156179, 17.649108886718864], [-101.918518066406179, 17.916076660156307], [-102.478088378906222, 17.97589111328125], [-103.500976562499957, 18.292297363281392], [-103.917480468749901, 18.748718261718892], [-104.992004394531165, 19.316284179687585], [-105.492980957031179, 19.946899414062614], [-105.731384277343665, 20.43408203125], [-105.397705078124915, 20.531677246093864], [-105.500610351562429, 20.81689453125], [-105.270690917968778, 21.0762939453125], [-105.265808105468665, 21.422119140625028], [-105.603088378906222, 21.87127685546875], [-105.693420410156151, 22.269104003906392], [-106.028686523437443, 22.773681640625028], [-106.909912109374972, 23.767883300781364], [-107.915405273437415, 24.548889160156307], [-108.401916503906207, 25.172302246093892], [-109.260192871093679, 25.5806884765625], [-109.444091796874957, 25.824890136718807], [-109.291625976562401, 26.442871093750057], [-109.801391601562415, 26.67608642578125], [-110.391723632812429, 27.162109375000028], [-110.640991210937401, 27.859924316406392], [-111.178894042968778, 27.941284179687557], [-111.759582519531179, 28.468078613281278], [-112.228210449218679, 28.954528808593864], [-112.271789550781165, 29.266906738281307], [-112.809509277343707, 30.021118164062614], [-113.163818359374972, 30.786926269531364], [-113.148681640624915, 31.171081542968892], [-113.871887207031278, 31.567687988281278], [-114.205688476562457, 31.524108886718864], [-114.776428222656179, 31.799682617187585], [-114.936706542968736, 31.393493652343892], [-114.771179199218665, 30.913696289062557], [-114.673889160156151, 30.1627197265625], [-114.330993652343693, 29.750488281250142], [-113.588806152343665, 29.061706542968835], [-113.424011230468665, 28.826293945312557], [-113.271911621093679, 28.754882812499972], [-113.140014648437415, 28.411315917968864], [-112.962280273437443, 28.425292968749972], [-112.761596679687415, 27.780273437500028], [-112.457885742187457, 27.525878906250057], [-112.244873046874901, 27.171875], [-111.616516113281165, 26.662902832031222], [-111.284606933593722, 25.732727050781278], [-110.987792968749943, 25.294677734375085], [-110.710021972656236, 24.826110839843892], [-110.655029296874929, 24.298706054687614], [-110.172790527343707, 24.265686035156364], [-109.771789550781222, 23.811279296875028], [-109.409118652343651, 23.364685058593892], [-109.433410644531207, 23.185729980468835], [-109.854187011718722, 22.818298339843778], [-110.031311035156165, 22.823120117187557], [-110.294982910156236, 23.431091308593722], [-110.949523925781179, 24.001098632812557], [-111.670593261718665, 24.484497070312642], [-112.182006835937472, 24.738525390625114], [-112.148986816406193, 25.470275878906278], [-112.300720214843707, 26.012084960937642], [-112.777282714843707, 26.32208251953125], [-113.464599609374972, 26.768310546875085], [-113.596679687499957, 26.639526367187557], [-113.848876953124915, 26.900085449218835], [-114.465698242187401, 27.142089843750028], [-115.055114746093693, 27.722717285156307], [-114.982177734374957, 27.798278808593835], [-114.570312499999901, 27.74151611328125], [-114.199279785156278, 28.115112304687614], [-114.161987304687528, 28.566101074218778], [-114.931823730468707, 29.279479980468807], [-115.518676757812429, 29.556274414062585], [-115.887390136718693, 30.180908203125057], [-116.258300781249972, 30.836486816406335], [-116.721496582031222, 31.635681152343864], [-117.127685546874901, 32.5352783203125], [-117.295898437499972, 33.046325683593835], [-117.943908691406207, 33.621276855468892], [-118.410583496093665, 33.74090576171875], [-118.519897460937528, 34.027893066406364], [-119.080993652343736, 34.078125], [-119.438781738281207, 34.348510742187614], [-120.367797851562401, 34.447082519531335], [-120.622802734374901, 34.608520507812642], [-120.744323730468665, 35.15692138671875], [-121.714599609374929, 36.161682128906392], [-122.547485351562457, 37.551879882812585], [-122.512023925781222, 37.783508300781222], [-122.953186035156222, 38.113708496093778], [-123.727111816406151, 38.951721191406364], [-123.865112304687401, 39.767089843750114], [-124.398010253906222, 40.313293457031307], [-124.178771972656151, 41.142089843750085], [-124.213684082031151, 41.999694824218892], [-124.532775878906151, 42.766113281250028], [-124.142089843750014, 43.708496093749972], [-123.898925781249957, 45.523498535156392], [-124.079589843749901, 46.864685058593835], [-124.395690917968693, 47.720275878906278], [-124.687194824218665, 48.184509277343778], [-124.566101074218778, 48.379699707031335], [-123.119995117187486, 48.04010009765625], [-122.587280273437401, 47.096130371093892], [-122.340026855468764, 47.360107421875142], [-122.5, 48.18011474609375], [-122.840026855468665, 49.000122070312557], [-122.974182128906151, 49.002685546875114], [-124.910217285156222, 49.984680175781278], [-125.624572753906179, 50.416687011718778], [-127.435607910156278, 50.830688476562614], [-127.992675781249957, 51.715881347656364], [-127.850280761718793, 52.329711914062472], [-129.129699707031165, 52.755493164062528], [-129.305175781249915, 53.561706542968778], [-130.514892578124886, 54.287719726562585], [-130.536071777343665, 54.802673339843807], [-131.085815429687443, 55.178894042968892], [-131.967224121093722, 55.497924804687585], [-132.249999999999886, 56.370117187500057], [-133.539184570312386, 57.178894042968864], [-134.078002929687386, 58.12310791015625], [-135.038208007812528, 58.187683105468835], [-136.627990722656278, 58.212280273437585], [-137.799987792968722, 58.500122070312614], [-139.8677978515625, 59.537902832031335], [-140.825195312500028, 59.727478027343835], [-142.574401855468665, 60.084472656250028], [-143.958801269531193, 59.999328613281392], [-145.925476074218693, 60.45867919921875], [-147.114379882812443, 60.884704589843864], [-148.224304199218693, 60.673095703125028], [-148.018005371093778, 59.978271484375057], [-148.570800781249886, 59.914306640625], [-149.727783203125057, 59.705688476562528], [-150.608215332031307, 59.368286132812557], [-151.716308593749886, 59.155883789062614], [-151.859375, 59.7451171875], [-151.409729003906193, 60.72589111328125], [-150.346923828124886, 61.033691406250142], [-150.621093750000057, 61.284484863281364], [-151.89581298828125, 60.727294921875057], [-152.578308105468693, 60.061706542968892], [-154.019104003906307, 59.350280761718722], [-153.287475585937472, 58.86468505859375], [-154.232482910156136, 58.146484375000085], [-155.307495117187443, 57.7279052734375], [-156.308288574218693, 57.422912597656222], [-156.556091308593807, 56.980102539062557], [-158.117187499999886, 56.463684082031392], [-158.433288574218693, 55.994079589843892], [-159.603271484375057, 55.566711425781278], [-160.2896728515625, 55.643676757812614], [-161.223022460937415, 55.364685058593722], [-162.237792968749972, 55.024291992187557], [-163.069396972656222, 54.68988037109375], [-164.785583496093665, 54.404296875000114], [-164.942199707031278, 54.572326660156364], [-163.848327636718693, 55.03948974609375], [-162.869995117187472, 55.348083496093722], [-161.80419921875, 55.895080566406364], [-160.563598632812557, 56.008117675781307], [-160.070495605468636, 56.418090820312642], [-158.684387207031193, 57.016723632812557], [-158.461120605468778, 57.2169189453125], [-157.722778320312443, 57.570129394531392], [-157.550292968749972, 58.328308105468892], [-157.041687011718693, 58.918884277343778], [-158.194702148437472, 58.615905761718864], [-158.517211914062443, 58.787902832031307], [-159.058593749999943, 58.424316406250028], [-159.711608886718807, 58.931518554687557], [-159.981201171874943, 58.57269287109375], [-160.35528564453125, 59.071105957031278], [-161.354980468750057, 58.6708984375], [-161.968811035156193, 58.671691894531278], [-162.054992675781136, 59.266906738281392], [-161.874084472656136, 59.63372802734375], [-162.518005371093636, 59.989685058593778], [-163.818298339843778, 59.798095703125114], [-164.66217041015625, 60.267517089843722], [-165.346374511718693, 60.50750732421875], [-165.350769042968693, 61.073913574218864], [-166.121398925781193, 61.500122070312528], [-165.734375, 62.075073242187614], [-164.919189453125, 62.633117675781307], [-164.5625, 63.146484374999972], [-163.753295898437386, 63.219482421875114], [-163.067199707031222, 63.059509277343722], [-162.260498046875057, 63.541870117187472], [-161.534423828124943, 63.45587158203125], [-160.772521972656193, 63.766113281250085], [-160.958312988281193, 64.222900390625085], [-161.51800537109375, 64.402893066406278], [-160.777709960937443, 64.788696289062472], [-161.391906738281136, 64.77728271484375], [-162.453002929687443, 64.559509277343778], [-162.757812500000057, 64.338684082031335], [-163.546386718749943, 64.55908203125], [-164.960815429687472, 64.447082519531392], [-166.425292968749972, 64.686706542968778], [-166.844970703125028, 65.088928222656222], [-168.110473632812415, 65.67010498046875], [-166.705200195312472, 66.088317871093864], [-164.474670410156165, 66.576721191406278], [-163.652526855468693, 66.576721191406278], [-163.78851318359375, 66.077270507812614], [-161.677795410156222, 66.116088867187528], [-162.489685058593636, 66.735473632812585], [-163.719726562499886, 67.116516113281307], [-164.430908203125, 67.616271972656335], [-165.390197753906165, 68.042907714843807], [-166.764404296874972, 68.358886718749972], [-166.204711914062472, 68.883117675781307], [-164.430786132812386, 68.915527343750057], [-163.168579101562472, 69.371093750000142], [-162.930480957031165, 69.85809326171875], [-161.908874511718807, 70.333312988281364], [-160.934814453124886, 70.447692871093807], [-159.0391845703125, 70.891723632812557], [-158.119689941406193, 70.824707031249972], [-156.580810546874886, 71.357910156250085], [-155.067810058593636, 71.147888183593807], [-154.34417724609375, 70.696472167968864], [-153.9000244140625, 70.890075683593864], [-152.210021972656136, 70.830078125], [-152.270019531249886, 70.600097656250057], [-150.739990234375057, 70.430114746093778], [-149.719970703125057, 70.530090332031307], [-147.613281249999943, 70.214111328125114], [-145.690002441406193, 70.120117187499972], [-144.919982910156193, 69.990112304687528], [-143.58941650390625, 70.152526855468892], [-142.072509765625028, 69.851928710937642], [-140.985900878906278, 69.712097167968864], [-139.120483398437415, 69.471130371093807], [-137.546386718749943, 68.990112304687557], [-136.503601074218665, 68.898071289062642], [-135.625671386718778, 69.315124511718807], [-134.414611816406307, 69.627502441406278], [-132.929199218749972, 69.505310058593835], [-131.431274414062386, 69.944519042968778], [-129.794677734374886, 70.193725585937642], [-129.107727050781222, 69.779296875], [-128.361511230468665, 70.012878417968864], [-128.138183593750028, 70.48388671875], [-127.447082519531207, 70.377319335937557], [-125.756286621093778, 69.480712890625085], [-124.424804687499972, 70.158508300781278], [-124.289611816406165, 69.399719238281222], [-123.061096191406179, 69.563720703125114], [-122.683410644531151, 69.855529785156364], [-121.472290039062528, 69.797912597656335], [-119.942810058593707, 69.377929687500114], [-117.602600097656222, 69.011291503906307], [-116.226379394531165, 68.841491699218778], [-115.246887207031278, 68.905883789062642], [-113.897888183593665, 68.398925781249972], [-115.304809570312472, 67.902709960937528], [-113.497192382812443, 67.68829345703125], [-110.797912597656193, 67.806091308593864], [-109.946105957031207, 67.981079101562528], [-108.880187988281207, 67.381530761718892], [-107.792419433593707, 67.887512207031335], [-108.812988281249901, 68.311706542968835], [-108.167175292968665, 68.653930664062642], [-106.950012207031179, 68.700073242187614], [-106.150024414062472, 68.800109863281278], [-105.342773437499915, 68.561279296875], [-104.337890624999929, 68.018127441406307], [-103.221069335937443, 68.097900390625114], [-101.454284667968778, 67.646911621093778], [-99.901977539062429, 67.805725097656392], [-98.443176269531222, 67.781677246093864], [-98.558593749999972, 68.4039306640625], [-97.669494628906165, 68.578674316406278], [-96.119873046874915, 68.239501953125], [-96.125793457031165, 67.293518066406307], [-95.489379882812472, 68.090698242187472], [-94.684997558593693, 68.063903808593807], [-94.232788085937443, 69.069091796875085], [-95.304077148437472, 69.685729980468807], [-96.471313476562415, 70.089904785156335], [-96.391113281249915, 71.194885253906364], [-95.208801269531222, 71.920471191406222], [-93.889892578124915, 71.760070800781222], [-92.878112792968722, 71.3187255859375], [-91.519592285156222, 70.1912841796875], [-92.406921386718636, 69.700073242187614], [-90.547119140624943, 69.497680664062585]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-114.16717, 73.12145], [-114.66634, 72.65277], [-112.44102, 72.9554], [-111.05039, 72.4504], [-109.92035, 72.961130000000111], [-109.00654, 72.63335], [-108.18835, 71.65089], [-107.68599, 72.06548], [-108.39639, 73.08953], [-107.51645, 73.23598], [-106.52259, 73.07601], [-105.40246, 72.67259], [-104.77484, 71.6984], [-104.46476, 70.99297], [-102.78537, 70.49776], [-100.98078, 70.024320000000102], [-101.08929, 69.58447], [-102.73116, 69.50402], [-102.09329, 69.11962], [-102.430239999999898, 68.75282], [-104.24, 68.91], [-105.96, 69.18], [-107.12254, 69.11922], [-109.0, 68.780000000000115], [-111.9668, 68.60446], [-113.3132, 68.53554], [-113.85496, 69.00744], [-115.22, 69.28], [-116.10794, 69.16821], [-117.34, 69.96], [-116.67473, 70.06655], [-115.13112, 70.2373], [-113.72141, 70.19237], [-112.4161, 70.366380000000106], [-114.35, 70.6], [-116.48684, 70.52045], [-117.9048, 70.54056], [-118.43238, 70.9092], [-116.11311, 71.30918], [-117.65568, 71.2952], [-119.40199, 71.55859], [-118.56267, 72.30785], [-117.86642, 72.70594], [-115.18909, 73.31459], [-114.16717, 73.12145]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-104.5, 73.42], [-105.38, 72.76], [-106.94, 73.46], [-106.6, 73.6], [-105.26, 73.64], [-104.5, 73.42]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-76.34, 73.102684989953019], [-76.25140380859375, 72.826385498046875], [-77.314437866210909, 72.855545043945284], [-78.39167022705081, 72.876655578613281], [-79.486251831054659, 72.742202758789091], [-79.775833129882841, 72.802902221679659], [-80.876098632812443, 73.333183288574304], [-80.833885192871065, 73.693183898925781], [-80.353057861328125, 73.759719848632784], [-78.064437866210938, 73.651931762695341], [-76.34, 73.102684989953019]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-86.56217851433405, 73.157447007938543], [-85.774371304044564, 72.534125881633827], [-84.850112474288153, 73.340278225387038], [-82.315590176101068, 73.750950832810673], [-80.600087653307554, 72.716543687624124], [-80.74894161652449, 72.061906643350682], [-78.770638597310779, 72.352173163534246], [-77.824623989559512, 72.749616604291049], [-75.605844692675731, 72.243678493937495], [-74.228616095664989, 71.767144273557903], [-74.099140794557712, 71.330840155717652], [-72.24222571479757, 71.556924546994509], [-71.200015428335121, 70.920012518997225], [-68.786054246684898, 70.525023708774341], [-67.914970465756852, 70.121947536897693], [-66.969033372654167, 69.186087348091888], [-68.805122850200547, 68.720198472764423], [-66.449866095633951, 68.067163397892017], [-64.862314419195144, 67.847538560651714], [-63.424934454996702, 66.928473212340663], [-61.851981370680505, 66.862120673277929], [-62.163176845942218, 66.160251369889608], [-63.918444383384099, 64.998668524832937], [-65.148860236253739, 65.426032619886684], [-66.721219041598545, 66.388041083432284], [-68.015016038673963, 66.262725735124491], [-68.141287400979252, 65.689789130304462], [-67.089646165623321, 65.108455105236999], [-65.732080451099677, 64.648405666758606], [-65.320167609301194, 64.38273712834615], [-64.669406297449598, 63.392926744227566], [-65.013803880458909, 62.674185085695996], [-66.275044725190384, 62.945098781986076], [-68.78318620469264, 63.745670071051819], [-67.369680752213156, 62.883965562584791], [-66.32829728866713, 62.280074774822054], [-66.165568203380161, 61.930897121825893], [-68.877366502544561, 62.330149237712902], [-71.023437059193924, 62.910708116295837], [-72.235378587519079, 63.397836005295176], [-71.886278449171215, 63.679989325608943], [-73.378306240518413, 64.193963121183913], [-74.834418911422688, 64.679075629323791], [-74.81850257027665, 64.389093329517948], [-77.709979824520133, 64.229542344816792], [-78.555948859354089, 64.572906399180141], [-77.897281053362036, 65.309192206474705], [-76.01827429879711, 65.326968899183242], [-73.959795294882639, 65.454764716240987], [-74.293883429649554, 65.811771348729309], [-73.944912482382648, 66.310578111426736], [-72.651167161739409, 67.284575507263867], [-72.926059943316005, 67.72692576768236], [-73.311617804645749, 68.069437160912827], [-74.843307257776729, 68.554627183701285], [-76.869100918266668, 68.894735622830353], [-76.228649054657268, 69.14776927354751], [-77.287369961237033, 69.769540106883284], [-78.168633999326516, 69.826487535268882], [-78.957242194316734, 70.166880194775416], [-79.492455003563663, 69.871807766388912], [-81.305470954091788, 69.743185126414346], [-84.944706183598583, 69.966634019644403], [-87.060003424817808, 70.26000112576537], [-88.681713223001424, 70.41074127876081], [-89.51341956252304, 70.762037665480904], [-88.467721116880767, 71.218185533321332], [-89.888151211287607, 71.222552191849957], [-90.205160285181933, 72.235074367960721], [-89.436576707705058, 73.129464219852451], [-88.408241543312812, 73.537888902471309], [-85.826151089201034, 73.803815823045227], [-86.56217851433405, 73.157447007938543]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-100.35642, 73.84389], [-99.16387, 73.63339], [-97.38, 73.76], [-97.12, 73.47], [-98.05359, 72.99052], [-96.54, 72.56], [-96.72, 71.66], [-98.35966, 71.27285], [-99.32286, 71.35639], [-100.01482, 71.73827], [-102.5, 72.51], [-102.48, 72.83], [-100.43836, 72.70588], [-101.54, 73.36], [-100.35642, 73.84389]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[143.60385, 73.21244], [142.08763, 73.20544], [140.038155, 73.31692], [139.86312, 73.36983], [140.81171, 73.765060000000119], [142.06207, 73.85758], [143.48283, 73.47525], [143.60385, 73.21244]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-93.196295539100134, 72.771992499473271], [-94.269046597047179, 72.024596259235977], [-95.409855516322665, 72.061880805134678], [-96.033745083382456, 72.940276801231903], [-96.018267991911017, 73.437429918095802], [-95.495793423223944, 73.862416897264183], [-94.503657599652428, 74.13490672473921], [-92.420012173211774, 74.100025132942278], [-90.509792853542507, 73.856732489712044], [-92.003965216829812, 72.96624420845859], [-93.196295539100134, 72.771992499473271]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-120.46, 71.4], [-123.09219, 70.901640000000128], [-123.62, 71.34], [-125.92896, 71.86868], [-125.59271, 72.19452], [-124.80729, 73.02256], [-123.94, 73.68], [-124.91775, 74.29275], [-121.53788, 74.44893], [-120.10978, 74.24135], [-117.55564, 74.18577], [-116.58442, 73.89607], [-115.51081, 73.47519], [-116.76794, 73.22292], [-119.22, 72.52], [-120.46, 71.820000000000135], [-120.46, 71.4]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[150.73167, 75.08406], [149.575925, 74.68892], [147.977465, 74.778355], [146.11919, 75.17298], [146.358485, 75.49682], [148.22223, 75.345845], [150.73167, 75.08406]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-93.612755906940407, 74.979997260224422], [-94.15690873897384, 74.592346503386864], [-95.608680589565694, 74.666863918751858], [-96.82093217648449, 74.92762319609659], [-96.28858740922982, 75.377828274223361], [-94.850819871789241, 75.647217515760985], [-93.977746548217937, 75.296489569796051], [-93.612755906940407, 74.979997260224422]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[145.086285, 75.562625], [144.3, 74.82], [140.61381, 74.84768], [138.95544, 74.61148], [136.97439, 75.26167], [137.51176, 75.94917], [138.831075, 76.13676], [141.471615, 76.09289], [145.086285, 75.562625]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-98.5, 76.72], [-97.735585, 76.256560000000121], [-97.704415, 75.74344], [-98.16, 75.0], [-99.80874, 74.89744], [-100.88366, 75.057360000000102], [-100.86292, 75.64075], [-102.50209, 75.5638], [-102.56552, 76.3366], [-101.48973, 76.30537], [-99.98349, 76.64634], [-98.57699, 76.58859], [-98.5, 76.72]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-108.21141, 76.20168], [-107.819429999999898, 75.84552], [-106.92893, 76.01282], [-105.881, 75.9694], [-105.70498, 75.47951], [-106.31347, 75.00527], [-109.7, 74.85], [-112.22307, 74.41696], [-113.74381, 74.39427], [-113.87135, 74.72029], [-111.79421, 75.1625], [-116.31221, 75.043430000000114], [-117.7104, 75.2222], [-116.34602, 76.19903], [-115.40487, 76.47887], [-112.59056, 76.14134], [-110.81422, 75.54919], [-109.0671, 75.47321], [-110.49726, 76.42982], [-109.5811, 76.79417], [-108.54859, 76.67832], [-108.21141, 76.20168]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[57.5356925799924, 70.720463975702245], [56.944979282463947, 70.632743231886678], [53.677375115784201, 70.762657782668555], [53.412016635965443, 71.206661688920292], [51.601894565645665, 71.474759019650406], [51.455753615124223, 72.014881089965144], [52.478275180883628, 72.229441636841045], [52.444168735570969, 72.774731350384769], [54.427613559797607, 73.627547512497671], [53.508289829325207, 73.749813951300155], [55.902458937407715, 74.627486477345428], [55.631932814359715, 75.081412258597169], [57.868643833248854, 75.609390367323215], [61.170044386647504, 76.251883450008222], [64.498368361270224, 76.439055487769281], [66.21097700385522, 76.809782213031212], [68.157059767534776, 76.939696763813004], [68.852211134725138, 76.544811306454534], [68.180572544227658, 76.233641669409025], [64.637326287703075, 75.737754625136233], [61.583507521414816, 75.260884507946884], [58.47708214705338, 74.309056301562833], [56.986785516188007, 73.333043524866156], [55.419335971910897, 72.37126760526607], [55.622837762276419, 71.540594794390415], [57.5356925799924, 70.720463975702245]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-94.684085862999467, 77.097878323058382], [-93.573921068073048, 76.776295884906091], [-91.605023159536643, 76.778517971494693], [-90.741845872749224, 76.449597479956736], [-90.969661424507905, 76.074013170059459], [-89.822237921899188, 75.847773749485725], [-89.187082892599904, 75.61016551380763], [-87.838276333349711, 75.566188869927231], [-86.379192267588593, 75.482421373182177], [-84.789625210290524, 75.699204006646511], [-82.753444586909978, 75.784315090631168], [-81.128530849924289, 75.713983466282002], [-80.057510952459154, 75.336848863415895], [-79.833932868148423, 74.923127346487206], [-80.457770758775922, 74.657303778777788], [-81.948842536125625, 74.44245901152442], [-83.228893602211429, 74.564027818490871], [-86.097452358733221, 74.410032050261151], [-88.150350307960338, 74.392307033985077], [-89.764722052758401, 74.51555532500123], [-92.422440965529432, 74.837757880341087], [-92.768285488642732, 75.386819973442158], [-92.88990597204176, 75.882655341282742], [-93.89382402217592, 76.31924367950063], [-95.962457445035824, 76.44138092722244], [-97.121378953829577, 76.7510777859477], [-96.745122850312271, 77.161388658345146], [-94.684085862999467, 77.097878323058382]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-116.198586595507351, 77.645286770326209], [-116.335813361458491, 76.876961575010654], [-117.106050584768795, 76.530031846819213], [-118.040412157038247, 76.48117178008718], [-119.899317586885701, 76.053213406061985], [-121.499995077126499, 75.900018622532798], [-122.854924486159078, 76.116542873835783], [-122.85492529360323, 76.116542873835783], [-121.157535360328254, 76.86450755482835], [-119.103938971821151, 77.512219957174636], [-117.570130784965968, 77.498318996888202], [-116.198586595507351, 77.645286770326209]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[106.970275878906278, 76.97430419921875], [107.240112304687585, 76.4801025390625], [108.153930664062585, 76.723327636718778], [111.077270507812528, 76.710083007812528], [113.331481933593778, 76.222290039062557], [114.134277343750028, 75.84771728515625], [113.885498046875085, 75.327880859375], [112.779296875000028, 75.031921386718722], [110.151306152343835, 74.476684570312557], [109.400085449218835, 74.180114746093722], [110.640075683593778, 74.040100097656222], [112.119323730468778, 73.787719726562528], [113.019470214843835, 73.976928710937528], [113.529724121093778, 73.3350830078125], [113.968872070312528, 73.594909667968807], [115.567871093750085, 73.752929687500028], [118.776306152343835, 73.58770751953125], [119.020080566406335, 73.1201171875], [123.200683593750085, 72.971313476562472], [123.257873535156278, 73.735107421875057], [125.380126953125085, 73.560119628906222], [126.976501464843835, 73.56549072265625], [128.591308593750028, 73.038696289062528], [129.051696777343835, 72.398681640625], [128.460083007812528, 71.9801025390625], [129.716125488281335, 71.193115234375], [131.288696289062528, 70.787109375000057], [132.253479003906278, 71.836303710937528], [133.857727050781335, 71.386474609375028], [135.562072753906278, 71.6552734375], [137.497680664062528, 71.34771728515625], [138.234130859375085, 71.62811279296875], [139.869873046875085, 71.4879150390625], [139.147888183593778, 72.41632080078125], [140.468078613281278, 72.849487304687528], [149.500122070312528, 72.200073242187557], [150.351318359375085, 71.606506347656307], [152.968872070312585, 70.84228515625], [157.006896972656278, 71.031494140625], [158.997924804687585, 70.86669921875], [159.830322265625028, 70.453308105468778], [159.708679199218835, 69.722106933593778], [160.940673828125028, 69.43731689453125], [162.279113769531335, 69.642089843750028], [164.052490234375085, 69.668273925781222], [165.940490722656335, 69.472106933593722], [167.835693359375085, 69.582702636718778], [169.577697753906335, 68.693908691406222], [170.816894531250085, 69.013671874999972], [170.008300781250028, 69.652893066406278], [170.453491210937528, 70.097106933593778], [173.643920898437642, 69.817504882812557], [175.724121093749972, 69.877319335937472], [178.600097656250028, 69.400085449218807], [180.000000000000142, 68.96372213254719], [180.000000000000142, 64.97958425748152], [179.992919921875142, 64.97430419921875], [178.707275390625028, 64.534912109375057], [177.411315917968892, 64.608276367187472], [178.313110351562642, 64.075927734374972], [178.908325195312528, 63.2520751953125], [179.370483398437585, 62.982727050781222], [179.486511230468892, 62.569091796875028], [179.228271484375142, 62.304077148437557], [177.364318847656278, 62.521911621093722], [174.569274902343778, 61.769287109375], [173.680114746093835, 61.652709960937528], [172.150085449218778, 60.950073242187472], [170.698486328125028, 60.3363037109375], [170.330871582031278, 59.88189697265625], [168.900512695312528, 60.573486328125], [166.295104980468778, 59.7886962890625], [165.840087890625028, 60.160095214843778], [164.876892089843778, 59.731689453125], [163.539306640625085, 59.86871337890625], [163.217102050781278, 59.211120605468807], [162.017272949218778, 58.2432861328125], [162.053100585937585, 57.839111328125], [163.191894531250085, 57.615112304687528], [163.057922363281278, 56.1593017578125], [162.129699707031335, 56.122314453125], [161.701477050781335, 55.285705566406307], [162.117492675781335, 54.855285644531278], [160.368896484375028, 54.344482421875028], [160.021728515625028, 53.202697753906222], [158.530883789062585, 52.95867919921875], [158.231323242187528, 51.942687988281307], [156.789916992187528, 51.011108398437528], [156.420104980468778, 51.700073242187557], [155.991882324218778, 53.158874511718778], [155.433715820312528, 55.381103515625], [155.914489746093835, 56.767883300781278], [156.758300781250085, 57.364685058593778], [156.810485839843835, 57.832092285156278], [158.364318847656278, 58.055725097656278], [160.150695800781278, 59.314880371093807], [161.872070312500028, 60.34307861328125], [163.669677734375085, 61.14093017578125], [164.473693847656278, 62.550720214843722], [163.258483886718778, 62.466308593749972], [162.657897949218778, 61.642517089843807], [160.121520996093835, 60.5443115234375], [159.302307128906278, 61.774108886718778], [156.720703125000028, 61.434509277343778], [154.218078613281278, 59.75830078125], [155.043884277343835, 59.14508056640625], [152.811889648437585, 58.883911132812557], [151.265686035156278, 58.780883789062528], [151.338073730468778, 59.504089355468807], [149.783691406250085, 59.655700683593778], [148.544921875000028, 59.164489746093807], [145.487304687500085, 59.33648681640625], [142.197875976562585, 59.040100097656278], [138.958496093750028, 57.088073730468722], [135.126281738281278, 54.729675292968778], [136.701721191406278, 54.603698730468778], [137.193481445312585, 53.977294921875], [138.164672851562528, 53.755126953125028], [138.804687500000028, 54.254699707031222], [139.901489257812528, 54.189697265625], [141.345275878906278, 53.0897216796875], [141.379272460937585, 52.238891601562528], [140.597473144531278, 51.23968505859375], [140.513122558593835, 50.045471191406222], [140.062072753906278, 48.446716308593807], [138.554687500000028, 46.999694824218807], [138.219726562500085, 46.307922363281222], [136.862304687500028, 45.143493652343722], [135.515319824218835, 43.98907470703125], [134.869506835937528, 43.398315429687528], [133.536926269531335, 42.811523437500028], [132.906311035156278, 42.79852294921875], [132.278076171875085, 43.284484863281222], [130.935913085937528, 42.55267333984375], [130.780090332031278, 42.220092773437472], [130.400085449218835, 42.28009033203125], [129.965881347656335, 41.941284179687528], [129.667480468750085, 41.60107421875], [129.705322265625085, 40.88287353515625], [129.188110351562528, 40.661926269531278], [129.010498046875085, 40.4854736328125], [128.633483886718835, 40.189880371093778], [127.967529296875085, 40.0255126953125], [127.533508300781278, 39.756896972656278], [127.502075195312585, 39.32391357421875], [127.385498046875028, 39.213500976562528], [127.783325195312528, 39.050903320312528], [128.349670410156335, 38.612304687500028], [129.212890625000028, 37.4324951171875], [129.460510253906335, 36.784301757812528], [129.468322753906278, 35.632080078125028], [129.091491699218778, 35.08251953125], [128.185913085937528, 34.890502929687528], [127.386474609375085, 34.4757080078125], [126.485717773437585, 34.39007568359375], [126.373901367187528, 34.934692382812472], [126.559326171875028, 35.6846923828125], [126.117492675781278, 36.72552490234375], [126.860290527343778, 36.893920898437472], [126.174682617187528, 37.749694824218778], [125.689086914062528, 37.94012451171875], [125.568481445312585, 37.752075195312528], [125.275329589843835, 37.669128417968778], [125.240112304687528, 37.857299804687557], [124.981079101562528, 37.948913574218778], [124.712280273437585, 38.108276367187528], [124.986083984375085, 38.548522949218807], [125.221923828125085, 38.665893554687528], [125.132873535156278, 38.84869384765625], [125.386718750000085, 39.387878417968722], [125.321105957031335, 39.551513671874972], [124.737487792968778, 39.6602783203125], [124.265686035156335, 39.92852783203125], [122.867675781250085, 39.637878417968778], [122.131530761718835, 39.17047119140625], [121.054687500000028, 38.89752197265625], [121.586120605468835, 39.36090087890625], [121.376892089843835, 39.75030517578125], [122.168701171875028, 40.422485351562528], [121.640502929687585, 40.94647216796875], [120.768676757812528, 40.593505859374972], [119.639709472656278, 39.898071289062557], [119.023498535156335, 39.252319335937528], [118.042724609375028, 39.204284667968722], [117.532714843750085, 38.7376708984375], [118.059692382812528, 38.0615234375], [118.878295898437528, 37.897277832031222], [118.911682128906278, 37.448486328124972], [119.702880859375028, 37.156494140625], [120.823486328125085, 37.870483398437472], [121.711303710937528, 37.481079101562557], [122.357910156250028, 37.454528808593778], [122.520080566406335, 36.930725097656222], [121.104309082031278, 36.651306152343722], [120.637084960937528, 36.111511230468807], [119.664672851562585, 35.60992431640625], [119.151306152343835, 34.909912109375], [120.227478027343778, 34.360473632812472], [120.620483398437528, 33.376708984375], [121.229125976562528, 32.460327148437528], [121.908081054687585, 31.69232177734375], [121.891906738281278, 30.949279785156222], [121.264282226562585, 30.676330566406307], [121.503479003906278, 30.14288330078125], [122.092102050781278, 29.832519531250057], [121.938476562500028, 29.01812744140625], [121.684509277343778, 28.22552490234375], [121.125671386718778, 28.13568115234375], [120.395507812500028, 27.053283691406278], [119.585510253906335, 25.74090576171875], [118.656921386718835, 24.5474853515625], [117.281677246093778, 23.62451171875], [115.890686035156335, 22.78289794921875], [114.763916015625085, 22.6680908203125], [114.152526855468835, 22.223876953125028], [113.806884765625028, 22.548278808593722], [113.241088867187585, 22.051513671875], [111.843688964843835, 21.55047607421875], [110.785522460937528, 21.397277832031222], [110.444091796875028, 20.34112548828125], [109.889892578125085, 20.282470703125], [109.627685546875028, 21.00830078125], [109.864501953125085, 21.395080566406222], [108.522888183593835, 21.71527099609375], [108.050292968750028, 21.552490234375028], [106.715087890625085, 20.696899414062472], [105.881713867187528, 19.752075195312557], [105.662109375000028, 19.058288574218807], [106.426879882812528, 18.004089355468807], [107.361877441406335, 16.697509765625], [108.269470214843778, 16.0797119140625], [108.877075195312585, 15.276672363281307], [109.335327148437585, 13.426086425781222], [109.200073242187585, 11.6668701171875], [108.366088867187585, 11.008300781250028], [107.220886230468778, 10.364501953124972], [106.405090332031278, 9.530883789062543], [105.158325195312585, 8.59967041015625], [104.795288085937528, 9.241088867187486], [105.076293945312585, 9.918518066406264], [104.334472656250028, 10.4866943359375], [103.497314453125085, 10.632690429687486], [103.090698242187528, 11.153686523437486], [102.585083007812528, 12.186706542968736], [101.687072753906278, 12.645874023437557], [100.831909179687528, 12.627075195312543], [100.978515625000028, 13.412719726562543], [100.097900390625028, 13.406921386718722], [100.018676757812585, 12.307128906250057], [99.478881835937528, 10.846496582031222], [99.153686523437585, 9.963073730468778], [99.222473144531335, 9.239318847656278], [99.873901367187585, 9.207885742187472], [100.279724121093835, 8.295288085937528], [100.459289550781335, 7.429687499999986], [101.017272949218835, 6.856872558593793], [101.623107910156278, 6.740722656250057], [102.141296386718835, 6.2216796875], [102.371276855468778, 6.128295898437486], [102.961730957031335, 5.524475097656278], [103.381286621093778, 4.855102539062543], [103.438720703125028, 4.181701660156278], [103.332092285156278, 3.726684570312514], [103.429504394531335, 3.38287353515625], [103.502502441406278, 2.791076660156293], [103.854675292968778, 2.5155029296875], [104.247924804687528, 1.631286621093778], [104.228881835937585, 1.2930908203125], [103.519714355468835, 1.226318359374986], [102.573730468750028, 1.967102050781293], [101.390686035156278, 2.760925292968778], [101.273681640625028, 3.270324707031236], [100.695495605468835, 3.939086914062528], [100.557495117187585, 4.767272949218793], [100.196716308593778, 5.3125], [100.306274414062528, 6.040710449218793], [100.085876464843778, 6.4644775390625], [99.690673828125028, 6.848327636718736], [99.519714355468835, 7.343505859375014], [98.988281250000028, 7.908081054687514], [98.503906250000028, 8.38232421875], [98.339721679687528, 7.794494628906264], [98.150085449218778, 8.350097656250014], [98.259277343750028, 8.973876953124986], [98.553527832031278, 9.933105468750014], [98.457275390625028, 10.675292968750043], [98.764526367187585, 11.441284179687543], [98.428283691406278, 12.0330810546875], [98.509704589843835, 13.122497558593778], [98.103698730468835, 13.6405029296875], [97.777709960937528, 14.8372802734375], [97.597106933593835, 16.1007080078125], [97.164672851562528, 16.928710937500028], [96.505920410156335, 16.427307128906222], [95.369323730468835, 15.714477539062528], [94.808471679687585, 15.80352783203125], [94.188903808593778, 16.038085937500057], [94.533508300781278, 17.277282714843778], [94.324890136718835, 18.213500976562528], [93.541076660156278, 19.366516113281307], [93.663330078125028, 19.727111816406222], [93.078308105468835, 19.855285644531222], [92.368530273437528, 20.6708984375], [92.082885742187585, 21.192321777343807], [92.025329589843778, 21.70172119140625], [91.834899902343778, 22.182922363281222], [91.417114257812585, 22.76507568359375], [90.496093750000028, 22.80511474609375], [90.587097167968778, 22.392883300781278], [90.272888183593835, 21.83648681640625], [89.847473144531335, 22.03912353515625], [89.702087402343778, 21.857116699218778], [89.418884277343835, 21.966308593750028], [89.032104492187528, 22.05572509765625], [88.888916015625085, 21.690673828125028], [88.208496093750085, 21.703308105468722], [86.975708007812528, 21.495483398437528], [87.033081054687585, 20.7432861328125], [86.499328613281278, 20.151672363281307], [85.060302734375085, 19.47869873046875], [83.941101074218778, 18.302124023437472], [83.189270019531278, 17.671325683593778], [82.192871093750028, 17.016723632812528], [82.191284179687585, 16.556701660156307], [81.692687988281335, 16.310302734375], [80.792114257812528, 15.952087402343722], [80.324890136718778, 15.8992919921875], [80.025085449218835, 15.136474609375028], [80.233276367187585, 13.835876464843722], [80.286315917968835, 13.006286621093736], [79.862487792968835, 12.056274414062514], [79.858093261718835, 10.357299804687514], [79.340515136718778, 10.308898925781307], [78.885498046875028, 9.546081542968736], [79.189697265625085, 9.216674804687486], [78.278076171875085, 8.933105468750028], [77.941284179687528, 8.252929687500014], [77.539916992187585, 7.965515136718793], [76.593078613281278, 8.899291992187486], [76.130126953125028, 10.2996826171875], [75.746520996093835, 11.30828857421875], [75.396118164062528, 11.78131103515625], [74.864929199218778, 12.741882324218778], [74.616699218750085, 13.992675781249986], [74.443908691406335, 14.617309570312557], [73.534301757812585, 15.99072265625], [73.119873046875028, 17.9287109375], [72.820922851562528, 19.208312988281222], [72.824523925781278, 20.41949462890625], [72.630676269531335, 21.356079101562472], [71.175292968750085, 20.757507324218778], [70.470520019531335, 20.877319335937557], [69.164123535156278, 22.08929443359375], [69.644897460937528, 22.450683593750028], [69.349670410156335, 22.843322753906307], [68.176696777343778, 23.692077636718778], [67.443725585937528, 23.944885253906278], [67.145507812500028, 24.663696289062557], [66.372924804687528, 25.425292968750057], [64.530517578125028, 25.237121582031278], [62.905700683593835, 25.218505859375028], [61.497497558593835, 25.07830810546875], [59.616088867187585, 25.380126953125], [58.525878906250028, 25.610107421875028], [57.397277832031335, 25.739929199218722], [56.970886230468835, 26.966125488281278], [56.492126464843835, 27.143310546875028], [55.723693847656278, 26.964721679687557], [54.715087890625028, 26.480712890625028], [53.493103027343835, 26.8125], [52.483703613281278, 27.580871582031278], [51.520874023437528, 27.865722656250028], [50.853088378906278, 28.81451416015625], [50.115112304687585, 30.147888183593807], [49.576904296875028, 29.9857177734375], [48.941284179687528, 30.317077636718807], [48.568115234375028, 29.9268798828125], [47.974487304687585, 29.97589111328125], [48.183288574218778, 29.534484863281307], [48.093872070312585, 29.306274414062557], [48.416076660156278, 28.552124023437472], [48.807678222656278, 27.689697265625057], [49.299682617187528, 27.461303710937557], [49.470886230468835, 27.110107421875], [50.152526855468778, 26.689697265624972], [50.213073730468778, 26.277099609375028], [50.113281250000085, 25.944091796875057], [50.239929199218835, 25.608093261718778], [50.527526855468778, 25.327880859375], [50.660705566406335, 24.999877929687472], [50.810119628906278, 24.754882812499972], [50.743896484375028, 25.482482910156278], [51.013488769531278, 26.007080078125], [51.286499023437585, 26.114685058593778], [51.589111328125085, 25.80108642578125], [51.606689453125028, 25.215698242187472], [51.389709472656335, 24.62750244140625], [51.579528808593778, 24.245483398437528], [51.757507324218835, 24.294128417968778], [51.794494628906278, 24.0198974609375], [52.577087402343778, 24.177490234375028], [53.404113769531278, 24.151306152343722], [54.008117675781335, 24.121887207031278], [54.693115234375085, 24.79791259765625], [55.439086914062585, 25.4390869140625], [56.070922851562528, 26.055480957031307], [56.362121582031278, 26.395874023437472], [56.485717773437528, 26.309082031249972], [56.391479492187528, 25.896118164062528], [56.261108398437528, 25.714721679687528], [56.396911621093778, 24.924682617187528], [56.845275878906335, 24.241699218750057], [57.403503417968778, 23.878723144531307], [58.137084960937528, 23.7479248046875], [58.729309082031278, 23.565673828125], [59.180480957031278, 22.992492675781222], [59.450073242187528, 22.660278320312528], [59.808105468750028, 22.53369140625], [59.806274414062528, 22.310485839843778], [59.442321777343835, 21.714477539062472], [59.282470703125085, 21.433898925781222], [58.861083984375028, 21.11407470703125], [58.488098144531335, 20.429077148437528], [58.034301757812585, 20.48150634765625], [57.826477050781278, 20.243103027343778], [57.665893554687528, 19.736083984374972], [57.788696289062528, 19.067687988281278], [57.694519042968835, 18.944702148437557], [57.234313964843778, 18.948120117187528], [56.609680175781278, 18.574279785156278], [56.512329101562528, 18.08709716796875], [56.283508300781278, 17.876098632812528], [55.661499023437585, 17.884277343750028], [55.270080566406335, 17.632324218750028], [55.274902343750028, 17.228271484375057], [54.791076660156335, 16.950683593749972], [54.239318847656278, 17.04510498046875], [53.570495605468778, 16.70770263671875], [53.108703613281335, 16.651123046874972], [52.385314941406335, 16.38250732421875], [52.191711425781335, 15.9384765625], [52.168273925781335, 15.597473144531278], [51.172485351562528, 15.175292968750028], [49.574707031250028, 14.708679199218722], [48.679321289062528, 14.003295898437514], [48.239074707031278, 13.948120117187557], [47.938903808593835, 14.007324218750043], [47.354492187500085, 13.592285156250043], [46.717102050781278, 13.399719238281264], [45.877685546875028, 13.347900390624986], [45.625122070312585, 13.29107666015625], [45.406494140625028, 13.026916503906264], [45.144287109375085, 12.953918457031222], [44.989685058593835, 12.69970703125], [44.494689941406278, 12.721679687500043], [44.175109863281278, 12.585876464843793], [43.483093261718778, 12.636901855468793], [43.222900390625085, 13.22088623046875], [43.251525878906335, 13.767700195312486], [43.088073730468778, 14.062683105468736], [42.892272949218835, 14.802307128906278], [42.604919433593835, 15.213317871093778], [42.805114746093835, 15.262084960937557], [42.702514648437528, 15.718872070312528], [42.823730468750028, 15.91168212890625], [42.779479980468778, 16.347900390625], [42.649719238281278, 16.774719238281222], [42.348083496093835, 17.075927734375], [42.270874023437528, 17.474670410156278], [41.754516601562585, 17.8331298828125], [41.221496582031278, 18.67169189453125], [40.939270019531278, 19.48651123046875], [40.247680664062528, 20.1746826171875], [39.801696777343835, 20.33892822265625], [39.139526367187528, 21.291870117187528], [39.023681640625085, 21.98687744140625], [39.066284179687528, 22.579711914062528], [38.492919921875028, 23.688476562499972], [38.023925781250085, 24.078674316406278], [37.483703613281335, 24.2855224609375], [37.154907226562528, 24.858520507812528], [37.209472656250028, 25.084472656249972], [36.931701660156335, 25.603088378906222], [36.639709472656278, 25.826293945312528], [36.249084472656278, 26.57012939453125], [35.640319824218835, 27.376525878906222], [35.130310058593835, 28.0634765625], [34.632324218750028, 28.058471679687528], [34.787902832031278, 28.607482910156278], [34.832275390625085, 28.957519531250057], [34.956115722656278, 29.356689453125028], [34.922729492187528, 29.501281738281278], [34.641723632812585, 29.099487304687528], [34.426696777343835, 28.344116210937557], [34.154479980468778, 27.823303222656278], [33.921508789062585, 27.648681640625028], [33.588073730468778, 27.971496582031278], [33.136901855468778, 28.417724609375028], [32.423278808593778, 29.851074218750028], [32.320495605468835, 29.760498046875028], [32.734924316406278, 28.705322265625], [33.348876953125028, 27.699890136718778], [34.104675292968778, 26.142272949218807], [34.473876953125085, 25.598693847656222], [34.795104980468835, 25.033874511718778], [35.692504882812585, 23.926696777343807], [35.493713378906278, 23.75250244140625], [35.526123046875028, 23.10247802734375], [36.690673828125028, 22.204895019531278], [36.866271972656335, 22.0001220703125], [37.188720703125085, 21.018920898437528], [36.969482421875085, 20.837524414062528], [37.114685058593778, 19.80810546875], [37.481872558593835, 18.614074707031222], [37.862670898437585, 18.367919921875057], [38.410095214843835, 17.998291015624972], [38.990722656250028, 16.840698242187557], [39.266113281250085, 15.9227294921875], [39.814270019531335, 15.435729980468778], [41.179321289062528, 14.4910888671875], [41.734924316406278, 13.921081542968764], [42.276916503906278, 13.344116210937514], [42.589721679687528, 13.00048828125], [43.081298828125085, 12.69970703125], [43.317871093750085, 12.390075683593778], [43.286499023437585, 11.97491455078125], [42.715881347656278, 11.735717773437486], [43.145324707031278, 11.462097167968736], [43.470703125000028, 11.277709960937528], [43.666687011718778, 10.864318847656307], [44.117919921875085, 10.445678710937486], [44.614318847656335, 10.442321777343722], [45.556884765625028, 10.698120117187486], [46.645507812500028, 10.816528320312543], [47.525695800781335, 11.127319335937528], [48.021728515625028, 11.193115234375028], [48.378906250000028, 11.375488281250057], [48.948303222656335, 11.410705566406236], [49.267883300781278, 11.430480957031278], [49.728698730468778, 11.578918457031222], [50.258911132812585, 11.679687500000028], [50.732116699218835, 12.021911621093736], [51.111328125000085, 12.024719238281278], [51.133911132812585, 11.748291015624972], [51.041503906250085, 11.166503906250028], [51.045288085937528, 10.640930175781264], [50.834289550781335, 10.279724121093722], [50.552490234375028, 9.198730468749972], [50.070922851562528, 8.08172607421875], [49.452697753906335, 6.804687500000028], [48.594482421875085, 5.339111328125057], [47.740905761718778, 4.219482421875043], [46.564880371093835, 2.855285644531236], [45.564086914062585, 2.045898437499986], [44.068298339843778, 1.052917480468807], [43.136108398437585, 0.29229736328125], [42.041687011718835, -0.919189453124943], [41.811096191406335, -1.446411132812472], [41.585083007812528, -1.6832275390625], [40.884887695312585, -2.082519531249986], [40.637878417968835, -2.499816894531207], [40.263122558593778, -2.573120117187514], [40.121276855468835, -3.277709960937528], [39.800109863281278, -3.68109130859375], [39.604919433593778, -4.346496582031278], [39.202270507812528, -4.676696777343722], [38.740478515625028, -5.90887451171875], [38.799682617187585, -6.475585937499957], [39.440124511718778, -6.840026855468736], [39.470092773437585, -7.099975585937472], [39.194702148437528, -7.703918457031207], [39.252075195312585, -8.0078125], [39.186523437500028, -8.485473632812472], [39.535888671875085, -9.112304687499957], [39.949707031250085, -10.098388671875028], [40.316711425781278, -10.317077636718693], [40.478515625000085, -10.765380859375014], [40.437316894531335, -11.761718749999957], [40.560913085937585, -12.639099121093707], [40.599670410156278, -14.201904296874986], [40.775512695312528, -14.691711425781236], [40.477294921875028, -15.406311035156236], [40.089294433593778, -16.100708007812514], [39.452697753906278, -16.720886230468764], [38.538330078125028, -17.101013183593764], [37.411071777343778, -17.586303710937486], [36.281311035156335, -18.659606933593764], [35.896484375000085, -18.842285156249972], [35.198486328125085, -19.55279541015625], [34.786499023437585, -19.783996582031278], [34.701904296875085, -20.497009277343722], [35.176086425781335, -21.254272460937514], [35.373474121093835, -21.840820312499972], [35.385925292968778, -22.140014648437486], [35.562683105468835, -22.090026855468722], [35.533874511718778, -23.070800781249986], [35.371887207031335, -23.535278320312528], [35.607482910156278, -23.706481933593707], [35.458679199218778, -24.12261962890625], [35.040710449218778, -24.478271484375028], [34.215881347656278, -24.816284179687472], [33.013305664062528, -25.357482910156207], [32.574707031250028, -25.727294921875028], [32.660278320312585, -26.148498535156236], [32.916076660156335, -26.215881347656207], [32.830078125000028, -26.742187500000014], [32.580322265625028, -27.470092773437486], [32.462280273437528, -28.301025390624986], [32.203491210937585, -28.752380371093707], [31.521118164062585, -29.25738525390625], [31.325683593750028, -29.401977539062514], [30.901672363281278, -29.909912109375014], [30.622924804687528, -30.423706054687472], [30.055725097656335, -31.140197753906222], [28.925476074218807, -32.171997070312486], [28.219726562500028, -32.771911621093793], [27.464721679687528, -33.226989746093778], [26.419494628906278, -33.614929199218764], [25.90972900390625, -33.6669921875], [25.78070068359375, -33.944580078124972], [25.17291259765625, -33.796875], [24.677917480468807, -33.987182617187457], [23.5941162109375, -33.79449462890625], [22.988281250000057, -33.9163818359375], [22.574279785156278, -33.864074707031207], [21.542907714843807, -34.258789062499957], [20.689086914062557, -34.417175292968778], [20.071289062500057, -34.795104980468764], [19.61651611328125, -34.819091796874972], [19.193298339843778, -34.462585449218793], [18.85528564453125, -34.444274902343793], [18.424682617187557, -33.997802734374986], [18.37750244140625, -34.136474609374986], [18.244506835937528, -33.867675781250028], [18.250122070312528, -33.281372070312457], [17.925292968750028, -32.611206054687528], [18.247924804687528, -32.429077148437457], [18.2216796875, -31.661621093750014], [17.56689453125, -30.725708007812543], [17.064514160156278, -29.878601074218736], [17.062927246093807, -29.875976562500028], [16.345092773437557, -28.576721191406222], [15.6019287109375, -27.821228027343722], [15.210510253906307, -27.090881347656207], [14.98968505859375, -26.117309570312486], [14.743286132812528, -25.392883300781222], [14.408081054687557, -23.853027343750014], [14.385681152343807, -22.656677246093707], [14.257690429687557, -22.111206054687486], [13.86871337890625, -21.698974609375028], [13.352478027343722, -20.872802734375014], [12.826904296875057, -19.673095703125028], [12.608703613281307, -19.045288085937528], [11.794921875000057, -18.069091796875014], [11.73431396484375, -17.301879882812528], [11.64007568359375, -16.673095703125014], [11.778686523437528, -15.793823242187528], [12.123718261718778, -14.878295898437528], [12.175720214843807, -14.449096679687486], [12.500122070312528, -13.547729492187528], [12.738525390625, -13.137878417968722], [13.31292724609375, -12.483581542968736], [13.63372802734375, -12.038574218749972], [13.738708496093807, -11.29779052734375], [13.686523437500057, -10.731079101562457], [13.387329101562528, -10.37359619140625], [13.121093750000057, -9.766906738281207], [12.875488281250028, -9.166870117187472], [12.9290771484375, -8.9591064453125], [13.236511230468722, -8.562622070312472], [12.933105468750028, -7.59649658203125], [12.728271484375057, -6.927124023437514], [12.227478027343778, -6.294372558593778], [12.322509765624972, -6.100097656250014], [12.182312011718807, -5.789916992187472], [11.915100097656222, -5.037902832031278], [11.093688964843807, -3.97882080078125], [10.066284179687557, -2.969482421875028], [9.405273437500028, -2.144287109374957], [8.798095703125057, -1.111328125000014], [8.830078125, -0.778991699218778], [9.048522949218807, -0.45928955078125], [9.291320800781278, 0.268676757812528], [9.492919921875028, 1.010070800781278], [9.305725097656307, 1.160888671874986], [9.6492919921875, 2.283874511718764], [9.7952880859375, 3.073486328124986], [9.40447998046875, 3.7344970703125], [8.948120117187557, 3.904113769531307], [8.744873046875057, 4.352294921875014], [8.488891601562528, 4.495727539062514], [8.50030517578125, 4.7720947265625], [7.462097167968807, 4.412109375000028], [7.082702636718778, 4.4647216796875], [6.698120117187472, 4.240722656250014], [5.898315429687528, 4.262512207031222], [5.3629150390625, 4.888122558593722], [5.03369140625, 5.611877441406236], [4.32568359375, 6.27069091796875], [3.57427978515625, 6.25830078125], [2.691711425781307, 6.258911132812543], [1.86529541015625, 6.142272949218778], [1.06011962890625, 5.92889404296875], [-0.50762939453125, 5.343505859374972], [-1.063598632812472, 5.000488281249986], [-1.964721679687472, 4.710510253906278], [-2.856079101562443, 4.994506835937514], [-3.311096191406193, 4.984313964843778], [-4.008789062499943, 5.179870605468736], [-4.64990234375, 5.168273925781278], [-5.834411621093693, 4.993713378906236], [-6.528686523437443, 4.705078125000028], [-7.518920898437472, 4.338317871093793], [-7.712097167968693, 4.36468505859375], [-7.974121093749972, 4.355895996093736], [-9.004821777343693, 4.832519531250057], [-9.913391113281222, 5.59368896484375], [-10.765380859374943, 6.140686035156222], [-11.438781738281193, 6.785888671875014], [-11.708190917968693, 6.860107421875043], [-12.428100585937528, 7.262878417968736], [-12.948974609375, 7.798706054687528], [-13.124023437499943, 8.163879394531293], [-13.24652099609375, 8.903076171875043], [-13.685180664062472, 9.494873046875], [-14.073974609374943, 9.886291503906278], [-14.330078125, 10.01568603515625], [-14.579711914062443, 10.214477539062543], [-14.693176269531278, 10.656311035156293], [-14.839477539062443, 10.876708984375043], [-15.13031005859375, 11.040527343750028], [-15.6641845703125, 11.45849609375], [-16.085205078124972, 11.524719238281293], [-16.314697265625, 11.806518554687528], [-16.308898925781193, 11.958679199218736], [-16.61376953125, 12.170898437500014], [-16.677429199218722, 12.3848876953125], [-16.84149169921875, 13.151489257812528], [-16.713684082031222, 13.595092773437486], [-17.126098632812528, 14.373474121093736], [-17.625, 14.729675292968722], [-17.185180664062472, 14.919494628906278], [-16.700683593749943, 15.621520996093807], [-16.463012695312443, 16.135070800781278], [-16.5496826171875, 16.673889160156278], [-16.270507812499972, 17.167114257812528], [-16.146301269531193, 18.108520507812557], [-16.25689697265625, 19.096679687500028], [-16.37762451171875, 19.593872070312557], [-16.27777099609375, 20.092529296875028], [-16.53631591796875, 20.567871093749972], [-17.06341552734375, 20.999877929687472], [-17.020385742187472, 21.42230224609375], [-16.97320556640625, 21.88568115234375], [-16.589111328124972, 22.1583251953125], [-16.261901855468693, 22.6793212890625], [-16.326416015625, 23.01788330078125], [-15.98260498046875, 23.7235107421875], [-15.426025390624972, 24.359130859375], [-15.089294433593722, 24.520324707031278], [-14.8245849609375, 25.103515625000028], [-14.800903320312472, 25.63629150390625], [-14.439880371093778, 26.254516601562557], [-13.7738037109375, 26.618896484375028], [-13.139892578125, 27.640075683593778], [-12.618774414062472, 28.038330078125028], [-11.68890380859375, 28.148681640625028], [-10.90087890625, 28.832275390625028], [-10.399597167968722, 29.09869384765625], [-9.564819335937443, 29.933715820312472], [-9.814697265624972, 31.177673339843807], [-9.434814453124943, 32.038085937499972], [-9.300720214843778, 32.564697265625057], [-8.657409667968722, 33.240295410156222], [-7.6541748046875, 33.697082519531307], [-6.912475585937472, 34.110473632812528], [-6.244323730468722, 35.145874023437528], [-5.929992675781222, 35.760070800781307], [-5.19378662109375, 35.755310058593807], [-4.591003417968693, 35.3306884765625], [-3.64007568359375, 35.39990234375], [-2.604309082031193, 35.179077148437528], [-2.169921874999943, 35.168518066406307], [-1.208618164062443, 35.71490478515625], [-0.127380371093722, 35.888671875], [0.503906250000057, 36.30133056640625], [1.466918945312557, 36.605712890624972], [3.161682128906307, 36.78387451171875], [4.815673828125, 36.865112304687472], [5.320129394531307, 36.716491699218778], [6.261901855468807, 37.1107177734375], [7.330505371093722, 37.1185302734375], [7.737121582031307, 36.885681152343778], [8.421081542968778, 36.946472167968722], [9.510070800781278, 37.35009765625], [10.210083007812528, 37.2301025390625], [10.180725097656307, 36.724121093750028], [11.028930664062557, 37.09210205078125], [11.10009765625, 36.900085449218778], [10.60009765625, 36.410095214843807], [10.59332275390625, 35.947509765625], [10.93951416015625, 35.6990966796875], [10.807922363281278, 34.833496093749972], [10.149719238281278, 34.330688476562528], [10.339721679687557, 33.785888671874972], [10.856872558593778, 33.768676757812528], [11.108520507812528, 33.293273925781278], [11.488891601562557, 33.137084960937557], [12.663330078125057, 32.792907714843807], [13.083312988281278, 32.878906250000028], [13.918701171875, 32.712097167968778], [15.2457275390625, 32.265075683593722], [15.713928222656278, 31.37628173828125], [16.6116943359375, 31.182312011718778], [18.0211181640625, 30.763488769531307], [19.086486816406278, 30.2664794921875], [19.574096679687528, 30.52587890625], [20.053283691406222, 30.985900878906222], [19.8203125, 31.751892089843722], [20.134094238281278, 32.23828125], [20.8544921875, 32.706909179687472], [21.543090820312557, 32.843322753906278], [22.895874023437557, 32.638488769531278], [23.236877441406278, 32.191528320312557], [23.609130859375057, 32.187316894531278], [23.927490234374972, 32.016723632812472], [24.921081542968778, 31.899475097656278], [25.164916992187528, 31.569274902343722], [26.495300292968778, 31.585693359375], [27.457702636718722, 31.321289062500057], [28.45050048828125, 31.02587890625], [28.91351318359375, 30.870117187500057], [29.683471679687585, 31.186889648437528], [30.095092773437585, 31.473510742187472], [30.976928710937528, 31.555908203124972], [31.688110351562585, 31.4296875], [31.960510253906278, 30.933715820312557], [32.192504882812585, 31.26031494140625], [32.993896484375085, 31.024108886718778], [33.773498535156335, 30.967529296875], [34.265502929687585, 31.219482421875], [34.556518554687528, 31.54888916015625], [34.488098144531335, 31.60552978515625], [34.752685546875028, 32.072875976562528], [34.955505371093835, 32.827514648437528], [35.098510742187528, 33.0806884765625], [35.126098632812585, 33.09088134765625], [35.482299804687585, 33.905517578125], [35.979675292968778, 34.610107421875028], [35.998474121093778, 34.6448974609375], [35.905090332031335, 35.410095214843722], [36.149902343750028, 35.821472167968778], [35.782104492187528, 36.275085449218722], [36.160888671875028, 36.650695800781278], [35.551086425781335, 36.56549072265625], [34.714477539062528, 36.795471191406278], [34.026916503906335, 36.220092773437528], [32.509277343750028, 36.107482910156278], [31.699707031250028, 36.644287109375028], [30.621704101562585, 36.677917480468722], [30.391113281250085, 36.263122558593778], [29.700073242187585, 36.144287109375028], [28.732910156250028, 36.6768798828125], [27.64129638671875, 36.65887451171875], [27.04888916015625, 37.653503417968807], [26.318298339843778, 38.208129882812528], [26.804687500000028, 38.985900878906222], [26.170898437500057, 39.463684082031222], [27.280090332031278, 40.420104980468807], [28.82012939453125, 40.460083007812472], [29.240112304687585, 41.2200927734375], [31.145874023437585, 41.087707519531222], [32.348083496093778, 41.736328125], [33.513305664062528, 42.019104003906278], [35.167724609375028, 42.040283203125028], [36.913085937500085, 41.335510253906278], [38.347717285156335, 40.948730468750057], [39.512695312500028, 41.102905273437528], [40.373474121093778, 41.013671875000057], [41.554077148437585, 41.535705566406278], [41.703308105468835, 41.963073730468722], [41.453491210937528, 42.64508056640625], [40.875488281250028, 43.013671875], [40.321472167968778, 43.128723144531278], [39.955078125000028, 43.435119628906278], [38.680114746093778, 44.280090332031278], [37.539123535156278, 44.657287597656307], [36.675476074218778, 45.24468994140625], [37.403320312500028, 45.404479980468778], [38.233093261718835, 46.24090576171875], [37.673706054687528, 46.636718750000028], [39.147705078125028, 47.044677734375028], [39.121276855468778, 47.263488769531307], [38.223693847656278, 47.102294921875028], [37.425109863281278, 47.022277832031278], [36.759887695312528, 46.698730468749972], [35.823730468750085, 46.645874023437528], [34.962280273437585, 46.2733154296875], [35.020874023437528, 45.651306152343807], [35.510070800781335, 45.410095214843778], [36.530090332031335, 45.470092773437528], [36.334716796875085, 45.11328125], [35.240112304687528, 44.940124511718778], [33.882507324218835, 44.36151123046875], [33.326477050781278, 44.564880371093807], [33.546875, 45.034912109375028], [32.454284667968835, 45.327514648437528], [32.630920410156335, 45.519287109375028], [33.588073730468778, 45.851684570312472], [33.298706054687528, 46.0806884765625], [31.744079589843835, 46.33349609375], [31.675292968750028, 46.706298828125], [30.748901367187528, 46.583129882812557], [30.377685546875085, 46.032470703125028], [29.603271484375085, 45.293273925781278], [29.626525878906278, 45.035522460937472], [29.1417236328125, 44.8203125], [28.837890625, 44.91387939453125], [28.558105468750028, 43.707519531250057], [28.039123535156278, 43.29327392578125], [27.673889160156307, 42.577880859375028], [27.996704101562557, 42.007507324218807], [28.115478515625, 41.6229248046875], [28.988525390625057, 41.2999267578125], [28.8065185546875, 41.054870605468807], [27.61907958984375, 40.9998779296875], [27.1925048828125, 40.690673828125057], [26.358093261718807, 40.152099609375028], [26.043273925781278, 40.617675781250028], [26.056884765625, 40.824096679687472], [25.447692871093778, 40.852478027343778], [24.9259033203125, 40.947082519531278], [23.714904785156307, 40.68707275390625], [24.408081054687528, 40.125122070312528], [23.900085449218778, 39.96209716796875], [23.343078613281222, 39.961120605468722], [22.8140869140625, 40.476074218750028], [22.626281738281307, 40.256530761718778], [22.84967041015625, 39.6593017578125], [23.350097656250028, 39.190124511718778], [22.97308349609375, 38.970886230468778], [23.530090332031222, 38.510070800781278], [24.025085449218778, 38.220092773437472], [24.040100097656307, 37.655090332031278], [23.115112304687528, 37.92010498046875], [23.410095214843807, 37.410095214843778], [22.77508544921875, 37.305114746093778], [23.154296875000028, 36.422485351562528], [22.490112304687472, 36.410095214843807], [21.670104980468778, 36.845092773437472], [21.295104980468778, 37.645080566406278], [21.120117187500028, 38.310302734375057], [20.730102539062472, 38.770080566406222], [20.21771240234375, 39.340270996093807], [20.15008544921875, 39.625122070312557], [19.980102539062528, 39.695129394531307], [19.960083007812557, 39.91510009765625], [19.4061279296875, 40.250915527343778], [19.319091796875028, 40.727294921875057], [19.403686523437528, 41.40948486328125], [19.540100097656307, 41.7200927734375], [19.37188720703125, 41.877685546875028], [19.1624755859375, 41.955078125], [18.882080078125, 42.281494140625057], [18.4500732421875, 42.4801025390625], [17.509887695312528, 42.850097656249972], [16.93011474609375, 43.210083007812557], [16.015502929687557, 43.507324218750028], [15.17449951171875, 44.243286132812528], [15.37628173828125, 44.317871093750057], [14.920288085937528, 44.738525390625028], [14.901672363281278, 45.076110839843778], [14.258728027343807, 45.233886718749972], [13.9522705078125, 44.802124023437528], [13.6571044921875, 45.137084960937557], [13.67950439453125, 45.484130859374972], [13.715087890625, 45.500305175781278], [13.93768310546875, 45.591125488281307], [13.1417236328125, 45.7366943359375], [12.328674316406278, 45.381896972656222], [12.383911132812557, 44.885498046875057], [12.261474609374972, 44.600524902343778], [12.589294433593778, 44.091491699218722], [13.52691650390625, 43.587707519531278], [14.0299072265625, 42.761108398437472], [15.142700195312557, 41.955078125], [15.926330566406307, 41.9613037109375], [16.169921874999972, 41.740295410156222], [15.8892822265625, 41.541076660156222], [16.785095214843778, 41.179687500000028], [17.519287109375, 40.877075195312528], [18.376708984375057, 40.355712890625028], [18.480285644531307, 40.168884277343778], [18.293518066406278, 39.810913085937472], [17.738525390624972, 40.277709960937528], [16.869689941406278, 40.442321777343778], [16.448730468750028, 39.795471191406222], [17.171508789062528, 39.424682617187472], [17.052917480468807, 38.902893066406278], [16.635070800781278, 38.843688964843807], [16.10107421875, 37.98590087890625], [15.684082031250028, 37.908874511718778], [15.688110351562557, 38.214721679687528], [15.892089843750057, 38.750915527343722], [16.109313964843778, 38.9644775390625], [15.718872070312528, 39.544128417968778], [15.413696289062528, 40.048278808593807], [14.998474121093778, 40.173095703125028], [14.703308105468778, 40.60467529296875], [14.060729980468778, 40.786499023437528], [13.62811279296875, 41.188293457031307], [12.888122558593807, 41.253112792968778], [12.106689453125057, 41.704528808593722], [11.191894531250028, 42.35552978515625], [10.5120849609375, 42.931518554687557], [10.200073242187528, 43.920104980468722], [9.7025146484375, 44.036315917968778], [8.888916015625, 44.36627197265625], [8.428710937500028, 44.2313232421875], [7.850891113281278, 43.76727294921875], [7.435302734375057, 43.693908691406222], [6.529296875000057, 43.128906250000028], [4.556884765625028, 43.39971923828125], [3.100524902343778, 43.075317382812528], [2.986083984375028, 42.473083496093807], [3.039489746093778, 41.892089843750028], [2.091918945312528, 41.226074218749972], [0.810485839843778, 41.014709472656278], [0.721313476562528, 40.678283691406222], [0.106689453125057, 40.124084472656307], [-0.278686523437443, 39.310119628906278], [0.111328125000028, 38.738525390625], [-0.46710205078125, 38.292480468749972], [-0.683410644531278, 37.642272949218722], [-1.438293457031222, 37.443115234375057], [-2.14642333984375, 36.674072265625057], [-3.415710449218693, 36.65887451171875], [-4.368896484374943, 36.677917480468722], [-4.995178222656193, 36.32470703125], [-5.377075195312472, 35.946899414062557], [-5.866394042968722, 36.0299072265625], [-6.236694335937472, 36.367675781249972], [-6.520202636718693, 36.94287109375], [-7.453674316406222, 37.097900390624972], [-7.8555908203125, 36.838317871093722], [-8.382812499999943, 36.978881835937472], [-8.898803710937472, 36.868896484375028], [-8.746093749999972, 37.651489257812528], [-8.840026855468693, 38.266296386718778], [-9.287475585937443, 38.358520507812528], [-9.526489257812472, 38.73748779296875], [-9.4468994140625, 39.39208984375], [-9.048278808593778, 39.755126953124972], [-8.977294921875, 40.159301757812472], [-8.768676757812443, 40.760681152343722], [-8.790771484375, 41.184326171875], [-8.990783691406193, 41.543518066406278], [-9.034790039062472, 41.880676269531222], [-8.984375, 42.5928955078125], [-9.392883300781222, 43.026672363281278], [-7.978210449218722, 43.74847412109375], [-6.754516601562443, 43.567871093750028], [-5.41180419921875, 43.574279785156278], [-4.3477783203125, 43.403503417968722], [-3.51751708984375, 43.455871582031222], [-1.901306152343722, 43.422912597656278], [-1.384216308593722, 44.022705078125057], [-1.19378662109375, 46.014892578125], [-2.225708007812443, 47.064514160156278], [-2.963195800781193, 47.5703125], [-4.491577148437472, 47.955078125000028], [-4.592285156249972, 48.684082031250057], [-3.295776367187472, 48.901672363281278], [-1.616516113281193, 48.644470214843778], [-1.933410644531193, 49.776489257812528], [-0.989379882812443, 49.347473144531222], [1.33868408203125, 50.1273193359375], [1.63909912109375, 50.94671630859375], [2.513488769531307, 51.14849853515625], [3.315124511718807, 51.345886230468807], [3.830322265625057, 51.6204833984375], [4.70611572265625, 53.0919189453125], [6.074279785156278, 53.510498046875028], [6.9052734375, 53.482299804687557], [7.100524902343778, 53.693908691406278], [7.936279296874972, 53.748291015624972], [8.1217041015625, 53.527893066406307], [8.80072021484375, 54.0208740234375], [8.572082519531307, 54.39569091796875], [8.526306152343807, 54.962890624999972], [8.120300292968778, 55.517700195312528], [8.090087890625057, 56.54010009765625], [8.256713867187472, 56.810119628906278], [8.543518066406278, 57.110107421875], [9.42449951171875, 57.172119140625], [9.775695800781307, 57.447875976562557], [10.580078125, 57.7301025390625], [10.546081542968807, 57.21588134765625], [10.250122070312528, 56.890075683593722], [10.370117187500028, 56.610107421875], [10.912292480468807, 56.45867919921875], [10.66790771484375, 56.08148193359375], [10.370117187500028, 56.19012451171875], [9.650085449218778, 55.4700927734375], [9.921875, 54.983093261718807], [9.939697265625028, 54.596679687499972], [10.950073242187557, 54.363708496093722], [10.93951416015625, 54.008728027343807], [11.956298828124972, 54.196472167968778], [12.518493652343722, 54.470520019531222], [13.647521972656307, 54.075500488281222], [14.119689941406307, 53.757080078124972], [14.802917480468807, 54.05072021484375], [16.363525390625, 54.513305664062528], [17.622924804687557, 54.851684570312557], [18.620910644531278, 54.682678222656278], [18.696289062500057, 54.438720703125028], [19.660705566406278, 54.42608642578125], [19.888488769531222, 54.866088867187557], [21.268493652343778, 55.190490722656278], [21.055908203125, 56.0311279296875], [21.09051513671875, 56.78387451171875], [21.581909179687557, 57.411926269531222], [22.52447509765625, 57.75347900390625], [23.3184814453125, 57.006286621093778], [24.1207275390625, 57.025695800781222], [24.312927246093778, 57.79351806640625], [24.429077148437528, 58.38348388671875], [24.061279296875057, 58.257507324218778], [23.426696777343807, 58.612670898437528], [23.339904785156307, 59.18731689453125], [24.604309082031222, 59.46588134765625], [25.86431884765625, 59.611083984375028], [26.949279785156307, 59.445922851562472], [27.981079101562472, 59.47552490234375], [29.11767578125, 60.028076171875], [28.070129394531307, 60.50347900390625], [26.255310058593722, 60.423889160156278], [24.496704101562528, 60.057312011718807], [22.86968994140625, 59.84649658203125], [22.2908935546875, 60.39190673828125], [21.322326660156307, 60.720275878906278], [21.544921874999972, 61.705322265625], [21.059326171875, 62.60748291015625], [21.53607177734375, 63.18988037109375], [22.442687988281307, 63.817871093749972], [24.730529785156278, 64.902282714843722], [25.3980712890625, 65.111511230468722], [25.29412841796875, 65.534484863281222], [23.903503417968778, 66.006896972656278], [22.183288574218807, 65.723876953124972], [21.213500976562557, 65.026123046875028], [21.369689941406278, 64.413696289062557], [19.7789306640625, 63.609680175781222], [17.847900390625057, 62.749511718750028], [17.11968994140625, 61.341308593750028], [17.831481933593807, 60.63671875], [18.787719726562528, 60.081909179687528], [17.869323730468778, 58.953918457031307], [16.829284667968807, 58.719909667968722], [16.447692871093778, 57.04107666015625], [15.879882812500057, 56.104309082031278], [14.666687011718778, 56.200927734375028], [14.100708007812557, 55.40789794921875], [12.942871093750057, 55.361877441406278], [12.6251220703125, 56.30712890625], [11.7880859375, 57.44189453125], [11.027282714843778, 58.856079101562472], [10.356689453125057, 59.46990966796875], [8.382080078125028, 58.31329345703125], [7.04888916015625, 58.078918457031307], [5.6658935546875, 58.588073730468778], [5.308288574218778, 59.663330078125028], [4.99212646484375, 61.971130371093778], [5.912902832031278, 62.614501953125057], [8.55352783203125, 63.454101562500028], [10.5277099609375, 64.486083984375057], [12.358276367187472, 65.879699707031222], [14.761291503906278, 67.81072998046875], [16.4359130859375, 68.563293457031278], [19.184082031250057, 69.817504882812557], [21.378479003906307, 70.255310058593778], [23.023681640625057, 70.202087402343778], [24.546691894531307, 71.030517578125], [26.370117187500057, 70.986328125000057], [28.16552734375, 71.18548583984375], [31.293518066406278, 70.453918457031222], [30.005493164062528, 70.186279296875028], [31.101074218750028, 69.558105468750057], [32.132690429687528, 69.905883789062528], [33.775512695312585, 69.301513671875], [36.514099121093835, 69.0634765625], [40.292480468750028, 67.932495117187472], [41.059875488281335, 67.457275390624972], [41.126098632812528, 66.791687011718807], [40.015930175781278, 66.266296386718778], [38.382873535156278, 65.99951171875], [33.918701171875028, 66.759704589843778], [33.184509277343835, 66.632507324218778], [34.814880371093778, 65.900085449218778], [34.943908691406278, 64.414489746093722], [36.231323242187585, 64.109497070312472], [37.012878417968835, 63.849914550781222], [37.142089843750085, 64.334716796875028], [36.517700195312585, 64.780273437500028], [37.176086425781278, 65.143310546875], [39.593505859375085, 64.520874023437528], [40.435729980468778, 64.764526367187528], [39.762695312500085, 65.496887207031307], [42.093078613281335, 66.476318359374972], [43.016113281250028, 66.418701171875057], [43.949890136718835, 66.069091796875057], [44.532287597656335, 66.756286621093778], [43.698486328125085, 67.35247802734375], [44.187927246093835, 67.950500488281307], [43.452880859375028, 68.570922851562528], [46.250122070312528, 68.250122070312528], [46.821289062500085, 67.689880371093722], [45.555297851562585, 67.566528320312528], [45.562072753906335, 67.010070800781278], [46.349121093750085, 66.667724609374972], [47.894287109375085, 66.884521484375], [48.138671875000028, 67.52252197265625], [50.227722167968835, 67.998718261718778], [53.717529296875085, 68.85748291015625], [54.471679687500085, 68.80828857421875], [53.485900878906278, 68.201293945312528], [54.726318359375085, 68.09710693359375], [55.442687988281278, 68.438720703125], [57.317077636718835, 68.466308593750028], [58.802124023437585, 68.880920410156222], [59.941528320312528, 68.27850341796875], [61.077880859375085, 68.940673828125028], [60.030090332031335, 69.520080566406222], [60.550109863281335, 69.850097656250028], [63.504089355468778, 69.547485351562528], [64.888122558593778, 69.234924316406222], [68.512084960937528, 68.092285156250028], [69.180725097656278, 68.61572265625], [68.164489746093778, 69.144287109375028], [68.135314941406278, 69.356506347656222], [66.930114746093835, 69.454711914062528], [67.259887695312585, 69.928710937499972], [66.724914550781278, 70.708923339843722], [66.694702148437528, 71.029113769531278], [68.540100097656335, 71.93450927734375], [69.196289062500028, 72.843505859375028], [69.940124511718778, 73.04010009765625], [72.587524414062528, 72.776306152343722], [72.796081542968778, 72.220092773437557], [71.848083496093835, 71.409118652343722], [72.470092773437585, 71.090270996093778], [72.791870117187585, 70.391113281250028], [72.564697265625085, 69.020874023437528], [73.667907714843835, 68.407897949218722], [73.238708496093778, 67.740478515625028], [71.280090332031335, 66.32012939453125], [72.423095703125028, 66.172729492187528], [72.820678710937585, 66.53271484375], [73.921081542968778, 66.789489746093807], [74.186523437500085, 67.2843017578125], [75.052124023437528, 67.760498046875028], [74.469299316406278, 68.329101562500057], [74.935913085937585, 68.989318847656278], [73.842285156250028, 69.071472167968722], [73.601928710937585, 69.627685546875], [74.399902343750028, 70.63189697265625], [73.101074218750028, 71.44708251953125], [74.890930175781278, 72.121276855468778], [74.659301757812528, 72.832275390624972], [75.158081054687528, 72.8551025390625], [75.683471679687528, 72.300476074218778], [75.289123535156278, 71.335693359375], [76.359130859375028, 71.15289306640625], [75.903076171875028, 71.874084472656278], [77.576721191406335, 72.267272949218778], [79.652099609375085, 72.320129394531278], [81.500122070312528, 71.750122070312557], [80.610717773437528, 72.582885742187472], [80.511108398437585, 73.648315429687557], [82.250122070312585, 73.850097656250028], [84.655273437500085, 73.805908203125], [86.822326660156335, 73.936889648437528], [86.009704589843835, 74.459716796874972], [87.166870117187585, 75.116516113281222], [88.315673828125085, 75.1439208984375], [90.260070800781278, 75.640075683593722], [92.900695800781278, 75.773315429687472], [93.234313964843778, 76.047302246093722], [95.860107421875085, 76.140075683593807], [96.678283691406278, 75.91552734375], [98.922485351562528, 76.446899414062472], [100.759704589843835, 76.430297851562472], [101.035278320312528, 76.86187744140625], [101.990905761718835, 77.2874755859375], [104.351684570312528, 77.697875976562528], [106.066711425781335, 77.3739013671875], [104.705078125000028, 77.127502441406278], [106.970275878906278, 76.97430419921875]], [[49.110290527343778, 41.28228759765625], [49.618896484375085, 40.572875976562528], [50.084899902343778, 40.52630615234375], [50.392883300781278, 40.256530761718778], [49.569274902343778, 40.176086425781222], [49.395324707031278, 39.399475097656278], [49.223327636718835, 39.049316406249972], [48.856506347656278, 38.81549072265625], [48.883300781250028, 38.320312500000057], [49.199707031250028, 37.5828857421875], [50.147888183593778, 37.37469482421875], [50.842285156250085, 36.872924804687557], [52.264099121093835, 36.700500488281307], [53.825927734375085, 36.965087890625], [53.921691894531335, 37.19891357421875], [53.735473632812528, 37.906127929687557], [53.880920410156278, 38.952087402343778], [53.101074218750085, 39.29071044921875], [53.357910156250028, 39.975280761718778], [52.694091796875028, 40.033691406249972], [52.915283203125085, 40.876525878906278], [53.858276367187585, 40.631103515625], [54.736877441406335, 40.951110839843807], [54.008300781250085, 41.551330566406278], [53.721679687500028, 42.123291015624972], [52.916687011718778, 41.86810302734375], [52.814697265625028, 41.135498046875], [52.502502441406278, 41.7833251953125], [52.446289062500028, 42.02728271484375], [52.692077636718778, 42.443908691406307], [52.501525878906278, 42.792297363281222], [51.342529296875028, 43.133117675781307], [50.891296386718835, 44.031127929687472], [50.339111328125085, 44.284118652343807], [50.305725097656335, 44.60992431640625], [51.278503417968835, 44.514892578125028], [51.316894531250028, 45.246093749999972], [52.167480468750028, 45.408508300781222], [53.040893554687585, 45.25909423828125], [53.220886230468778, 46.234680175781222], [53.042724609375028, 46.853088378906278], [52.042114257812585, 46.804687499999972], [51.192077636718835, 47.048706054687528], [50.034118652343835, 46.609130859375], [49.101318359375028, 46.399475097656222], [48.645507812500085, 45.806274414062557], [47.675903320312585, 45.641479492187557], [46.682128906250028, 44.609313964843807], [47.590881347656278, 43.660278320312528], [47.492492675781335, 42.9866943359375], [48.584472656250085, 41.80889892578125], [49.110290527343778, 41.28228759765625]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-93.840003017944014, 77.519997260234504], [-94.295608283245173, 77.491342678528781], [-96.169654100309998, 77.555111395976979], [-96.436304490936124, 77.83462921824372], [-94.422577277386296, 77.820004787905077], [-93.720656297565966, 77.634331366680328], [-93.840003017944014, 77.519997260234504]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-110.186938035913073, 77.697014879050386], [-112.051191169058484, 77.409228827616943], [-113.534278937619163, 77.732206529441243], [-112.724586758253949, 78.051050116682035], [-111.264443325630864, 78.152956041161644], [-109.854451870547109, 77.996324774884926], [-110.186938035913073, 77.697014879050386]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[24.72412, 77.85385], [22.49032, 77.44493], [20.72601, 77.67704], [21.41611, 77.93504], [20.8119, 78.25463], [22.88426, 78.45494], [23.28134, 78.07954], [24.72412, 77.85385]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-109.6631457182026, 78.60197256134569], [-110.881314256618964, 78.406919867659923], [-112.542091437615184, 78.407901719873507], [-112.525890876091694, 78.550554511215154], [-111.500010342233409, 78.849993598130567], [-110.963660651476033, 78.804440823065306], [-109.6631457182026, 78.60197256134569]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-95.830294969449255, 78.056941229963343], [-97.309842902398032, 77.850597235821795], [-98.124289313534092, 78.082856960757681], [-98.552867804746569, 78.458105373845029], [-98.631984422585447, 78.871930243638388], [-97.337231411512533, 78.831984361476685], [-96.754398769908704, 78.765812689927088], [-95.559277920294505, 78.418314520980289], [-95.830294969449255, 78.056941229963343]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-100.060191820052239, 78.32475434031582], [-99.67093909381353, 77.907544664207506], [-101.303940192452927, 78.018984890444898], [-102.949808722733053, 78.343228664860305], [-105.176132778731557, 78.38033234324584], [-104.210429450277175, 78.677420152491777], [-105.419580451258554, 78.918335679836531], [-105.492289191493256, 79.30159393992912], [-103.529282396237846, 79.165349026191649], [-100.825158047268729, 78.800461737778704], [-100.060191820052239, 78.32475434031582]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.5
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[105.07547, 78.30689], [99.43814, 77.921], [101.2649, 79.23399], [102.08635, 79.34641], [102.837815, 79.28129], [105.37243, 78.713340000000102], [105.07547, 78.30689]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[18.25183, 79.70175], [21.54383, 78.95611], [19.02737, 78.5626], [18.47172, 77.82669], [17.59441, 77.63796], [17.1182, 76.80941], [15.91315, 76.77045], [13.76259, 77.38035], [14.669560000000104, 77.73565], [13.170600000000121, 78.02493], [11.22231, 78.8693], [10.4445300000001, 79.65239], [13.17077, 80.01046], [13.71852, 79.66039], [15.14282, 79.67431], [15.52255, 80.01608], [16.99085, 80.05086], [18.25183, 79.70175]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[25.44762535981198, 80.407340399894593], [27.407505730913499, 80.056405748200547], [25.924650506298178, 79.517833970854554], [23.024465773213706, 79.400011705229105], [20.075188429451885, 79.566823228667261], [19.897266473070999, 79.842361965647513], [18.462263624757895, 79.859880276194417], [17.368015170977543, 80.318896186027104], [20.455992059010697, 80.59815562613224], [21.907944777115489, 80.357679348462057], [22.919252557067438, 80.657144273593502], [25.44762535981198, 80.407340399894593]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 1,
      "min_zoom": 1.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[51.136186557831394, 80.54728017854103], [49.79368452332082, 80.415427761548216], [48.894411248577597, 80.339566758943789], [48.754936557821765, 80.175468248200929], [47.586119012244211, 80.010181179515342], [46.502825962109711, 80.247246812654367], [47.072455275262968, 80.559424140129465], [44.846958042181228, 80.589809882317098], [46.799138624871233, 80.771917629713727], [48.318477410684665, 80.784009914869955], [48.522806023966751, 80.514568996900238], [49.097189568890855, 80.753985907708341], [50.039767693894674, 80.918885403151819], [51.52293297710375, 80.69972565380192], [51.136186557831394, 80.54728017854103]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[99.93976, 78.88094], [97.75794, 78.7562], [94.97259, 79.04474500000012], [93.31288, 79.426500000000118], [92.5454, 80.14379], [91.181070000000119, 80.34146], [93.77766, 81.0246], [95.940895, 81.2504], [97.88385, 80.746975000000106], [100.186655, 79.780135], [99.93976, 78.88094]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-87.02, 79.66], [-85.81435, 79.3369], [-87.18756, 79.0393], [-89.03535, 78.28723], [-90.80436, 78.21533], [-92.87669, 78.34333], [-93.95116, 78.75099], [-93.935739999999868, 79.11373], [-93.14524, 79.380100000000112], [-94.974, 79.37248], [-96.07614, 79.705020000000104], [-96.70972, 80.157770000000113], [-96.01644, 80.60233], [-95.32345, 80.907290000000103], [-94.29843, 80.97727], [-94.73542, 81.20646], [-92.40984, 81.2573900000001], [-91.13289, 80.72345], [-89.45, 80.509322033898286], [-87.81, 80.320000000000107], [-87.02, 79.66]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-68.5, 83.106321516765831], [-65.82735, 83.028010000000137], [-63.68, 82.9], [-61.85, 82.628600000000148], [-61.89388, 82.36165], [-64.334, 81.92775], [-66.75342, 81.72527], [-67.65755, 81.50141], [-65.48031, 81.50657], [-67.84, 80.9], [-69.4697, 80.616830000000135], [-71.18, 79.8], [-73.2428, 79.63415], [-73.88, 79.430162204802087], [-76.90773, 79.32309], [-75.52924, 79.19766], [-76.22046, 79.01907], [-75.393449999999888, 78.52581], [-76.343539999999876, 78.18296], [-77.88851, 77.89991], [-78.36269, 77.50859], [-79.75951, 77.20968], [-79.61965, 76.98336], [-77.91089, 77.022045], [-77.88911, 76.777955], [-80.56125, 76.17812], [-83.17439, 76.45403], [-86.11184, 76.29901], [-87.6, 76.42], [-89.49068, 76.47239], [-89.6161, 76.95213], [-87.76739, 77.178330000000102], [-88.26, 77.9], [-87.65, 77.970222222222304], [-84.97634, 77.53873], [-86.34, 78.18], [-87.96192, 78.37181], [-87.15198, 78.75867], [-85.37868, 78.9969], [-85.09495, 79.34543], [-86.50734, 79.73624], [-86.93179, 80.25145], [-84.19844, 80.20836], [-83.408695652173833, 80.100000000000108], [-81.84823, 80.46442], [-84.1, 80.58], [-87.59895, 80.51627], [-89.36663, 80.85569], [-90.2, 81.26], [-91.36786, 81.5531], [-91.58702, 81.89429], [-90.1, 82.085], [-88.93227, 82.11751], [-86.97024, 82.27961], [-85.5, 82.652273458057039], [-84.260005, 82.6], [-83.18, 82.32], [-82.42, 82.86], [-81.1, 83.02], [-79.30664, 83.13056], [-76.25, 83.172058823529397], [-75.71878, 83.06404], [-72.83153, 83.23324], [-70.665765, 83.169780758382927], [-68.5, 83.106321516765831]]]
    }
  }, {
    "type": "Feature",
    "properties": {
      "featurecla": "Land",
      "scalerank": 0,
      "min_zoom": 0.0
    },
    "geometry": {
      "type": "Polygon",
      "coordinates": [[[-27.10046, 83.51966], [-20.84539, 82.72669], [-22.69182, 82.34165], [-26.51753, 82.297650000000118], [-31.9, 82.200000000000102], [-31.39646, 82.02154], [-27.85666, 82.131780000000106], [-24.84448, 81.78697000000011], [-22.903279999999882, 82.0931700000001], [-22.07175, 81.734490000000108], [-23.169609999999892, 81.15271], [-20.62363, 81.52462], [-15.768179999999887, 81.91245], [-12.77018, 81.71885], [-12.20855, 81.29154], [-16.28533, 80.58004], [-16.85, 80.35], [-20.04624, 80.17708], [-17.73035, 80.12912], [-18.9, 79.4], [-19.70499, 78.75128], [-19.67353, 77.63859], [-18.47285, 76.98565], [-20.03503, 76.94434], [-21.67944, 76.62795], [-19.83407, 76.09808], [-19.59896, 75.24838], [-20.66818, 75.15585], [-19.37281, 74.29561], [-21.59422, 74.223820000000103], [-20.43454, 73.81713], [-20.76234, 73.46436], [-22.17221, 73.30955], [-23.56593, 73.30663], [-22.31311, 72.62928], [-22.29954, 72.18409], [-24.27834, 72.59788], [-24.79296, 72.330200000000104], [-23.44296, 72.08016], [-22.13281, 71.46898000000013], [-21.753559999999879, 70.66369], [-23.53603, 70.471], [-24.30702, 70.85649], [-25.54341, 71.43094], [-25.20135, 70.75226], [-26.36276, 70.22646], [-23.72742, 70.18401], [-22.34902, 70.12946], [-25.02927, 69.2588], [-27.74737, 68.47046], [-30.67371, 68.125030000000123], [-31.77665, 68.12078000000011], [-32.81105, 67.73547], [-34.20196, 66.67974], [-36.35284, 65.9789], [-37.04378, 65.93768], [-38.37505, 65.69213], [-39.81222, 65.45848], [-40.66899, 64.83997], [-40.68281, 64.13902], [-41.1887, 63.48246], [-42.81938, 62.68233], [-42.41666, 61.900930000000102], [-42.86619, 61.07404], [-43.3784, 60.09772], [-44.7875, 60.03676], [-46.26364, 60.85328], [-48.26294, 60.85843], [-49.23308, 61.40681], [-49.90039, 62.38336], [-51.63325, 63.62691], [-52.14014, 64.278420000000125], [-52.27659, 65.1767], [-53.66166, 66.09957], [-53.30161, 66.8365], [-53.969109999999887, 67.18899], [-52.9804, 68.35759], [-51.475359999999881, 68.729580000000112], [-51.08041, 69.14781], [-50.87122, 69.9291], [-52.013585, 69.574925], [-52.55792, 69.42616], [-53.45629, 69.283625], [-54.68336, 69.61003], [-54.75001, 70.28932], [-54.35884, 70.821315], [-53.431315, 70.835755], [-51.39014, 70.569780000000122], [-53.10937, 71.20485], [-54.00422, 71.54719], [-55.0, 71.406536967272501], [-55.83468, 71.65444], [-54.71819, 72.58625], [-55.326339999999874, 72.95861], [-56.12003, 73.649770000000103], [-57.32363, 74.71026], [-58.59679, 75.09861], [-58.585159999999888, 75.51727000000011], [-61.26861, 76.10238], [-63.391649999999885, 76.1752], [-66.06427, 76.13486], [-68.50438, 76.06141], [-69.66485, 76.37975], [-71.40257, 77.008570000000105], [-68.77671, 77.32312], [-66.76397, 77.37595], [-71.04293, 77.63595], [-73.297, 78.04419], [-73.159379999999885, 78.43271], [-69.37345, 78.91388], [-65.7107, 79.39436], [-65.3239, 79.75814], [-68.02298, 80.11721], [-67.15129, 80.51582], [-63.68925, 81.2139600000001], [-62.23444, 81.3211], [-62.65116, 81.77042], [-60.28249, 82.033630000000102], [-57.20744, 82.19074], [-54.13442, 82.19962], [-53.04328, 81.88833], [-50.39061, 82.43883], [-48.00386, 82.06481], [-46.59984, 81.9859450000001], [-44.523, 81.6607], [-46.9007, 82.19979], [-46.76379, 82.62796], [-43.40644, 83.22516], [-39.89753, 83.18018], [-38.62214, 83.54905], [-35.08787, 83.64513], [-27.10046, 83.51966]]]
    }
  }]
};
},{}],"index.ts":[function(require,module,exports) {
"use strict";

var __importStar = this && this.__importStar || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) {
    if (Object.hasOwnProperty.call(mod, k)) result[k] = mod[k];
  }
  result["default"] = mod;
  return result;
};

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

exports.__esModule = true;

var geo = __importStar(require("d3-geo"));

var d3_selection_1 = require("d3-selection"); // import world from './assets/world-110m.json'


var ne_110m_land_json_1 = __importDefault(require("./assets/ne_110m_land.json")); // const svg = document.getElementById('map')


var app = d3_selection_1.select('#app');
app.append('h1').attr('class', 'title').text('Coordinates of Site and Sound');
var svg = app.append('div').attr('id', 'map').style('width', 400).append('svg');
svg.attr('viewBox', '0 0 960 600').style('width', '100%').style('height', 'auto');
svg.append("path").datum({
  type: "FeatureCollection",
  features: ne_110m_land_json_1["default"].features
}).attr("d", geo.geoPath().projection(geo.geoEquirectangular())) // .attr("d", geo.geoPath().projection(geo.geoMercator()) as any)
// .attr("d", geo.geoPath().projection(geo.geoNaturalEarth1()) as any)
// .attr("d", geo.geoPath().projection(geo.geoEqualEarth()) as any)
.style('fill', 'white').style('stroke', '#000').style('stroke-width', '0.5px'); // svg.selectAll("path")
//   .data(land)
//   .enter().append("path")
//     .attr("d", geo.geoPath());
// svg.append('path')
//   .datum(land)
//   .enter().append('path')
//     .attr('d', geo.geoPath().projection(geo.geoMercator()))
// svg.selectAll("path")
//   // .data(topojson.feature(world, world.objects.land).features)
//   // .datum(topojson.mesh(world, world.objects.land))
//   .enter().append("path")
//   .attr("d", path);
},{"d3-geo":"../node_modules/d3-geo/src/index.js","d3-selection":"../node_modules/d3-selection/src/index.js","./assets/ne_110m_land.json":"assets/ne_110m_land.json"}],"../node_modules/parcel-bundler/src/builtins/hmr-runtime.js":[function(require,module,exports) {
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
  var ws = new WebSocket(protocol + '://' + hostname + ':' + "62232" + '/');

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
},{}]},{},["../node_modules/parcel-bundler/src/builtins/hmr-runtime.js","index.ts"], null)
//# sourceMappingURL=/src.77de5100.js.map