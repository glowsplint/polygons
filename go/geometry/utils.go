package geometry

import (
	dec "github.com/shopspring/decimal"
)

var zero = dec.Zero
var one = dec.NewFromInt(1)
var two = dec.NewFromInt(2)
var half = dec.NewFromFloat(0.5)

var zeroZero = NewPoint(zero, zero)
var zeroOne = NewPoint(zero, one)
var oneZero = NewPoint(one, zero)
var oneOne = NewPoint(one, one)
var halfHalf = NewPoint(half, half)

var Start = oneZero
var End = NewPoint(one.Neg(), zero)
