"""Read-level mean quality (the ONT ``qs`` value) as a filter on aligned reads.

Dorado writes the mean read quality in the ``qs`` tag, computed as
``-10 * log10(mean(10 ** (-q / 10)))`` over the base qualities, and nagini's
pre-alignment read filter (nanostream) uses the same definition.  Applying it
again on the aligned BAM lets a stricter cut-off be tried without aligning
again.  It can only be stricter than the pre-alignment filter, because reads
below that one are already gone.

The tag is used when the BAM has it.  Otherwise the value is recomputed from the
base qualities the way dorado computes ``qs``: from base 60 of the read as
sequenced (the adapter and barcode region is left out), or over every base for
a read split from another (``pi`` tag).  That reproduces ``qs`` for 99 % of
reads and matches nanostream's own fallback; averaging every base instead read
MinKNOW's Q20 reads as Q16.  A
supplementary alignment that is hard-clipped carries only part of the read's
qualities, so its value cannot be recomputed; such an alignment is kept and
counted as unknown rather than judged on a fragment, unless another alignment of
the same read was already found below the cut-off, in which case the whole read
goes.

Standard library only.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Set, Tuple

QS_TAG = "qs"
# Dorado leaves the first bases of a read out of its qs.
BASECALLER_QSCORE_START = 60


def mean_read_quality(qualities: Sequence[int]) -> Optional[float]:
    """-10*log10 of the mean error probability, or None without qualities."""

    if qualities is None or len(qualities) == 0:
        return None
    mean_error = sum(10.0 ** (-q / 10.0) for q in qualities) / len(qualities)
    return -10.0 * math.log10(mean_error) if mean_error > 0 else float(max(qualities))


def basecaller_mean_read_quality(qualities: Sequence[int], split_read: bool = False) -> Optional[float]:
    """The read's qs as dorado computes it, from qualities in sequenced order."""

    if qualities is None or len(qualities) == 0:
        return None
    if split_read or len(qualities) <= BASECALLER_QSCORE_START:
        return mean_read_quality(qualities)
    return mean_read_quality(list(qualities)[BASECALLER_QSCORE_START:])


class ReadQualityFilter:
    """Keeps reads whose mean quality is at least ``min_qs``.  0 switches it off."""

    def __init__(self, min_qs: float = 0.0):
        self.min_qs = float(min_qs or 0.0)
        self._computed: Dict[Tuple, Optional[float]] = {}
        self.dropped: Set[str] = set()  # names of reads below the cut-off
        self.unknown: Set[Tuple] = set()

    @property
    def enabled(self) -> bool:
        return self.min_qs > 0

    @staticmethod
    def _key(aln) -> Tuple:
        return (aln.query_name, bool(aln.is_supplementary), aln.reference_start, bool(aln.is_read1))

    def value(self, aln) -> Optional[float]:
        """The read's mean quality, or None when it cannot be known."""

        if aln.has_tag(QS_TAG):
            try:
                return float(aln.get_tag(QS_TAG))
            except (TypeError, ValueError):
                pass
        key = self._key(aln)
        if key in self._computed:
            return self._computed[key]
        qualities = aln.query_qualities
        complete = qualities is not None and len(qualities) == aln.infer_read_length()
        result = None
        if complete:
            # An alignment stores a reverse-strand read's qualities reversed.
            ordered = list(qualities)[::-1] if getattr(aln, "is_reverse", False) else list(qualities)
            result = basecaller_mean_read_quality(ordered, split_read=aln.has_tag("pi"))
        self._computed[key] = result
        return result

    def passes(self, aln) -> bool:
        if not self.enabled:
            return True
        if aln.query_name in self.dropped:
            return False
        value = self.value(aln)
        if value is None:
            self.unknown.add(self._key(aln))
            return True
        if value < self.min_qs:
            self.dropped.add(aln.query_name)
            return False
        return True

    def summary(self) -> str:
        if not self.enabled:
            return "read quality filter off"
        return (
            f"min_read_qs={self.min_qs:g}: {len(self.dropped)} reads dropped, "
            f"{len(self.unknown)} kept without a known qs"
        )
