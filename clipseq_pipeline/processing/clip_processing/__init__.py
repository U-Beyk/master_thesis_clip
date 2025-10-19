'''
clip_processing
===============

Package for processing and filtering CLIP data.

This package exposes:

- "ClipEntry" from "clip_entry"
- "ClipProcessor" from "clip_processing"

author: U.B.
'''

from .clip_entry import ClipEntry
from .clip_processor import ClipProcessor

__all__ = ["ClipEntry", "ClipProcessor"]