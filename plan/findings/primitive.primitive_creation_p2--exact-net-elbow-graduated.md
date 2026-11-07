# Exact-net small elbow: graduated seeded cases

Addendum to the [original failure record](primitive.primitive_creation_p2--exact-net-elbow-unites-into-an-invalid-body-below-radius-0-1.md).

The follow-up in #2840, fixing #3087, graduates the original generated cases
`primitive.primitive_creation_p2/gen0092` and `gen0189`. Their actual seeded
inputs use radius 0.001, both legs 0.01, and construction tolerance 1e-7.
At checkpoint `ed6119602`, both bodies pass Level2 validity with no failed checks.
The measured volume is 6.283172176379765e-8, compared with
`pi * 0.001^2 * 0.02 = 6.283185307179586e-8`; relative error
2.089831698288911e-6 is within the unchanged 1e-4 volume assertion.

The actual runner regenerated both correct golden records and their known-failure
entries were removed in `1ffa9a25e`. These are generated golden cases, with no
existing kernel-case gallery entry. The original failure record remains unchanged
as historical evidence.

The [consolidated endpoint finding](curve.curve_p1--projected-endpoints-small-mitres.md#original-seeded-elbow-cases)
contains the exact inputs, reproducible replay helper, measured results, and
matching diagnostic renders before and after the endpoint fix. Later projected
selection review work and its validation are recorded there separately; this
addendum reports the explicit seeded-body checkpoint above.
