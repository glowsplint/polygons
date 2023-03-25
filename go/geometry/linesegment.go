package geometry

import (
	"errors"

	dec "github.com/shopspring/decimal"
)

type LineSegment struct {
	P1 *Point `json:"p1"`
	P2 *Point `json:"p2"`
}

func NewLineSegment(a, b *Point) *LineSegment {
	if a.X.Equal(b.X) && a.Y.Equal(b.Y) {
		return nil
	}

	if a.X.LessThan(b.X) {
		return &LineSegment{a, b}
	}

	if a.X.Equal(b.X) && a.Y.LessThan(b.Y) {
		return &LineSegment{a, b}
	}

	return &LineSegment{b, a}
}

// SquaredL2 returns the Squared L2 distance of a line segment
func (ls *LineSegment) SquaredL2() (dec.Decimal, error) {
	if ls == nil {
		return zero, errors.New("line segment is nil")
	}
	a := ls.P1.X.Sub(ls.P2.X).Pow(two)
	b := ls.P1.Y.Sub(ls.P2.Y).Pow(two)
	return a.Add(b), nil
}

// Direction returns true if P2 is closer to the End than P1
func (ls *LineSegment) Direction() (bool, error) {
	p1ToEnd, err := NewLineSegment(ls.P1, End).SquaredL2()
	if err != nil {
		return false, err
	}
	p2ToEnd, err := NewLineSegment(ls.P2, End).SquaredL2()
	if err != nil {
		return false, err
	}
	direction := p1ToEnd.GreaterThan(p2ToEnd)
	if p1ToEnd.Equal(p2ToEnd) {
		return false, errors.New("both points are equidistant")
	}
	return direction, nil
}
