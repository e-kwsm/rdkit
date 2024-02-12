#!/usr/bin/env python3
import unittest

import rdmolfiles  # FIXME


class Test(unittest.TestCase):
    def testCML(self):
        m = rdmolfiles.MolFromCMLBlock(
            """<?xml version="1.0"?>
<cml>
  <molecule>
    <atomArray>
      <atom id="a0" elementType="Du"/>
    </atomArray>
    <!--
    <bondArray>
    </bondArray>
    -->
  </molecule>
</cml>
            """
        )


def main():
    unittest.main()


if __name__ == "__main__":
    main()
