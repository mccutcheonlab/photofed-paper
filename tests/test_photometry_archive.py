"""Tests for restoring figure input files."""

import hashlib
import io
import json
import re
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from download_data import SEPARATE_FILES, download_file, download_inputs, extract_archive
from photometry_files import ARCHIVE_NAME, PICKLE_NAMES
class PhotometryArchiveTests(unittest.TestCase):
    """Exercise restoration and validation with synthetic input bytes."""

    def setUp(self) -> None:
        """Create isolated input and output directories."""
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.archive = self.root / ARCHIVE_NAME
        self.destination = self.root / "restored"

    def write_archive(self) -> None:
        """Create an archive with one synthetic member for every input."""
        with zipfile.ZipFile(self.archive, "w") as archive:
            for name in PICKLE_NAMES:
                archive.writestr(f"data/{name}", name.encode())

    def test_round_trip_and_existing_files_are_preserved(self) -> None:
        """An archive restores all selected inputs but never overwrites them."""
        self.write_archive()
        extract_archive(self.archive, self.destination)
        for name in PICKLE_NAMES:
            self.assertEqual((self.destination / name).read_bytes(), name.encode())
        with self.assertRaises(FileExistsError):
            extract_archive(self.archive, self.destination)
        self.assertEqual((self.destination / PICKLE_NAMES[0]).read_bytes(), PICKLE_NAMES[0].encode())

    def test_rejects_unexpected_archive_members(self) -> None:
        """Reject extra and path-traversal members before extracting anything."""
        with zipfile.ZipFile(self.archive, "w") as archive:
            for name in PICKLE_NAMES:
                archive.writestr(f"data/{name}", b"synthetic")
            archive.writestr("../unexpected", b"synthetic")
        with self.assertRaises(ValueError):
            extract_archive(self.archive, self.destination)
        self.assertFalse(self.destination.exists())

    def test_download_verifies_zenodo_checksum(self) -> None:
        """Use Zenodo metadata to reject corrupted archive bytes."""
        payload = b"synthetic archive"
        checksum = hashlib.md5(payload, usedforsecurity=False).hexdigest()
        metadata = {
            "files": [{
                "key": ARCHIVE_NAME,
                "checksum": f"md5:{checksum}",
                "links": {"self": "https://zenodo.org/archive"},
            }]
        }
        responses = [io.BytesIO(json.dumps(metadata).encode()), io.BytesIO(payload)]
        with patch("download_data.urllib.request.urlopen", side_effect=responses):
            download_file("123", self.archive)
        self.assertEqual(self.archive.read_bytes(), payload)

        metadata["files"][0]["checksum"] = "md5:" + "0" * 32
        responses = [io.BytesIO(json.dumps(metadata).encode()), io.BytesIO(payload)]
        with patch("download_data.urllib.request.urlopen", side_effect=responses):
            with self.assertRaisesRegex(ValueError, "checksum"):
                download_file("123", self.root / "bad.zip")

    def test_downloads_and_installs_all_inputs(self) -> None:
        """Restore the pickles and separate files from verified downloads."""
        self.write_archive()
        payloads = {ARCHIVE_NAME: self.archive.read_bytes()}
        payloads.update({name: name.encode() for name in SEPARATE_FILES})
        entries = [
            {
                "key": name,
                "checksum": "md5:" + hashlib.md5(content, usedforsecurity=False).hexdigest(),
                "links": {"self": f"https://zenodo.org/files/{name}"},
            }
            for name, content in payloads.items()
        ]

        def fake_open(url: str) -> io.BytesIO:
            if url == "https://zenodo.org/api/records/123":
                return io.BytesIO(json.dumps({"files": entries}).encode())
            return io.BytesIO(payloads[url.rsplit("/", 1)[-1]])

        with patch("download_data.urllib.request.urlopen", side_effect=fake_open):
            download_inputs("123", self.destination)

        for name in (*PICKLE_NAMES, *SEPARATE_FILES):
            self.assertEqual((self.destination / name).read_bytes(), name.encode())
        with patch("download_data.urllib.request.urlopen") as request:
            with self.assertRaises(FileExistsError):
                download_inputs("123", self.destination)
            request.assert_not_called()

    def test_failed_separate_download_leaves_no_inputs(self) -> None:
        """A missing or invalid separate file must not install partial inputs."""
        self.write_archive()
        payload = self.archive.read_bytes()
        entries = [{
            "key": ARCHIVE_NAME,
            "checksum": "md5:" + hashlib.md5(payload, usedforsecurity=False).hexdigest(),
            "links": {"self": "https://zenodo.org/files/archive"},
        }]

        def fake_open(url: str) -> io.BytesIO:
            if url == "https://zenodo.org/api/records/123":
                return io.BytesIO(json.dumps({"files": entries}).encode())
            return io.BytesIO(payload)

        with patch("download_data.urllib.request.urlopen", side_effect=fake_open):
            with self.assertRaisesRegex(ValueError, SEPARATE_FILES[0]):
                download_inputs("123", self.destination)
        self.assertFalse(self.destination.exists())

    def test_separate_files_match_notebook_inputs(self) -> None:
        """The renamed files and notebook input paths stay consistent."""
        notebook_dir = Path(__file__).resolve().parents[1] / "notebooks"
        required = {
            "fig2_food_and_bw.ipynb": SEPARATE_FILES[:1],
            "fig3_pellet_behaviour.ipynb": SEPARATE_FILES[1:],
            "figSx_retrieval_latencies.ipynb": SEPARATE_FILES[2:],
        }
        for notebook_name, names in required.items():
            notebook = json.loads((notebook_dir / notebook_name).read_text(encoding="utf-8"))
            text = "\n".join(
                "".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"
            )
            for name in names:
                self.assertIn(f'"{name}"', text)
            self.assertIn("DATA_DIR", text)
            if notebook_name != "fig2_food_and_bw.ipynb":
                self.assertIn("DATAFOLDER = DATA_DIR", text)

    def test_manifest_matches_figure_inputs(self) -> None:
        """Keep the photometry manifest aligned with figure notebook inputs."""
        notebook_dir = Path(__file__).resolve().parents[1] / "notebooks"
        sources = [
            "fig4_fr1_photo.ipynb",
            "fig5_progratio_photo.ipynb",
            "figSx_fr1_photo_pellet_offset.ipynb",
            "figSx_progratio_photo_offset.ipynb",
        ]
        notebooks = [json.loads((notebook_dir / name).read_text(encoding="utf-8")) for name in sources]
        texts = [
            "\n".join("".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code")
            for notebook in notebooks
        ]
        referenced = {name for text in texts for name in re.findall(r'"([^"/\\]+\.pickle)"', text)}
        self.assertEqual(referenced, set(PICKLE_NAMES))
        self.assertEqual(len(PICKLE_NAMES), len(set(PICKLE_NAMES)))
        self.assertIn(
            '"W3_3_NR_Sucrose_FR1_ttl_pellet_offset.pickle", '
            '"pellettype": "sucrose", "dietphase": "NR2"',
            texts[2],
        )


if __name__ == "__main__":
    unittest.main()
