# Projected face votes must preserve certificate refusals

Issue: #3107. Construction: `topology::brep::p10::faces_to_delete`.

The existing raw-face fixtures have one face wholly below the diagonal cutter
segment from (-1,0,0) to (1,1,0). The reference point (-1.5,-1,0) is also on its
right side. Thus KeepPoint retains the face and DeletePoint removes it, or the
operation refuses when it cannot certify a required projection.

The rectangle has x in [-.05,.05], y in [.39,.41]; its distance below the
line y=(x+1)/2 is at least .065. The lens's two quadratic boundary curves have
control points below that line except for their endpoints on it. Bernstein
convexity proves their interiors lie below it. These are geometric checks,
independent of a golden's previous output.

## Measured failure

At knot gaps 1.5e-12 and 1.9e-12, the rectangle's vertex and midpoint side
queries return NotWithinTolerance. The lens's endpoint queries return ordinary
on-curve ambiguity, while its midpoint queries return NotWithinTolerance.
Nevertheless both whole-face operations returned success with KeepPoint deleting
the face and DeletePoint retaining it. The four gap-1e-8 controls return the
geometrically correct KeepPoint/ DeletePoint face counts 1 / 0.

The vertex/midpoint errors and final trim outputs were replayed dynamically;
the precise final ray witness has not been instrumented. The source independently
shows all three vote sites used `ignored`, which erased certificate failures
before the ray ambiguity rule could choose the opposite reference side.

## Construction change

All three vote sites preserve the typed NotWithinTolerance error, including its
original payload. Other statuses retain the existing policy: ordinary on-curve
ambiguity remains an absent vote; cancellation, preparation and unported errors
still propagate. The reference vote already propagates errors directly.
No geometry, knot, solver tolerance or global error policy changes.

`projected_side_classification` replays the eight narrow-gap combinations and
four wider-gap controls through the public ProjectAndTrim API. The empty imprint
subset isolates classification: on a refused vote the original face remains.
This does not promise rollback of earlier imprint edits for a nonempty subset.

## Validation

The public regression failed against the parent construction before this caller
change: a narrow-gap operation returned success where the test required refusal.
The linked parent included the final-selection and directed-span minimum fixes;
its library SHA256 was
`b9697942d8f25fc6e4cf1d6f0b4d7203b9b926266573d6290cdb3ae859502ab9`.

Commands for the final source:

```sh
cargo test -p solidcore --test projected_side_classification
cargo test -p solidcore --lib projected_vote_preserves_refusal_and_control_signals
```

Green validation is pending. No golden has been blessed.
