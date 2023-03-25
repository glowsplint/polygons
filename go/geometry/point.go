package geometry

import (
	dec "github.com/shopspring/decimal"
)

type Point struct {
	X dec.Decimal `json:"x"`
	Y dec.Decimal `json:"y"`
}

func NewPoint(x, y dec.Decimal) *Point {
	return &Point{x, y}
}

// Equal returns true if both points are nil, or their coordinates are equal
func (a *Point) Equal(b *Point) bool {
	if a == nil && b == nil {
		return true
	}
	if a == nil || b == nil {
		return false
	}
	return a.X.Equal(b.X) && a.Y.Equal(b.Y)
}

// NewPointFromIntersection returns a point from the intersection of two line segments
// Constructor for Point via the intersection of two line segments using Cramer's rule
//
// References:
// https://stackoverflow.com/a/20679579
// https://observablehq.com/@toja/line-box-intersection
func NewPointFromIntersection(L1, L2 *LineSegment) *Point {
	a1 := L1.P2.X.Sub(L1.P1.X)
	b1 := L2.P1.X.Sub(L2.P2.X)
	c1 := L2.P1.X.Sub(L1.P1.X)
	a2 := L1.P2.Y.Sub(L1.P1.Y)
	b2 := L2.P1.Y.Sub(L2.P2.Y)
	c2 := L2.P1.Y.Sub(L1.P1.Y)

	d := a1.Mul(b2).Sub(b1.Mul(a2))
	dX := c1.Mul(b2).Sub(b1.Mul(c2))
	dY := a1.Mul(c2).Sub(c1.Mul(a2))

	if d.Equal(zero) {
		return nil
	}

	s := dX.Div(d)
	t := dY.Div(d)

	if !(isInBetween(zero, s, one) && isInBetween(zero, t, one)) {
		return nil
	}

	x := (one.Sub(s)).Mul(L1.P1.X).Add(s.Mul(L1.P2.X))
	y := (one.Sub(s)).Mul(L1.P1.Y).Add(s.Mul(L1.P2.Y))
	return NewPoint(x, y)
}

func isInBetween(left, x, right dec.Decimal) bool {
	return left.LessThanOrEqual(x) && x.LessThanOrEqual(right)
}
