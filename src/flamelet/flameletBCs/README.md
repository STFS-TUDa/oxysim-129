# Flamelet Boundary Conditions
↑ [OxySim-129](./../../../README.md)

[boundaryLookUpFvPatchScalarField.H](./boundaryLookUp/boundaryLookUpFvPatchScalarField.H)
[boundaryLookUpFvPatchScalarField.C](./boundaryLookUp/boundaryLookUpFvPatchScalarField.C)

## `boundaryLookup`

A `fixedValue` BC (usually on `ha`) that just tags a patch for the non-adiabatic boundary lookup described in [`flameletThermo`](./../flameletThermo#non-adiabatic-boundary-lookup) — it stores a `lookupType` and otherwise behaves like `fixedValue`. It doesn't compute its own value: `flameletThermo::setBoundaryLookupType()` scans for patches of this type at start-up, and `flameletThermo::boundaryLookUp()` writes the looked-up `ha` into them on every `update()`.

```
wall
{
    type       boundaryLookup;
    lookupType "specialBoundaryLookup";  // or "standardBoundaryLookup" (default)
    value      uniform -254452.7113;
}
```

`lookupType` must match across every `boundaryLookup` patch in a case (`flameletThermo` enforces this).
