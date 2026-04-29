import os
import argparse
import json
import ROOT

root_file = "/scratchnvme/store/mc/RunIII2024Summer24NanoAODv15/DY2Mu_2Jets_MLL-105to160-April2026_IreneFakes_Oversampling9_FlashSim/260421_142939/0000/tree_10.root"

spec = {"samples": {"sample": {
         "trees": ["Events"],
         "files": ["/scratchnvme/store/mc/RunIII2024Summer24NanoAODv15/DY2Mu_2Jets_MLL-105to160-April2026_IreneFakes_Oversampling9_FlashSim/260421_142939/0000/tree_10.root"]
      }
   }
}
if not os.path.exists("spec.json"):
   with open("spec.json", "w") as f:
      json.dump(spec, f, indent=4)

parser = argparse.ArgumentParser()
parser.add_argument("-j", type=int, default=1, help="Number of threads to use (default: 1)")

args = parser.parse_args()
if args.j > 1:
   print(f"Starting test with implicit multi-threading enabled. Using {args.j} threads.")
   ROOT.EnableImplicitMT(args.j)

rdf = ROOT.RDF.Experimental.FromSpec("spec.json")
rdf = (
   rdf.DefinePerSample(
      "range_key",
      "std::to_string(static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().first)) + ':' + "
      "std::to_string(static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().second))",
   )
)

raw_keys = rdf.Take["std::string"]("range_key").GetValue()

ranges = []
for key in raw_keys:
   left, right = str(key).split(":", 1)
   ranges.append((int(left), int(right)))

# print("Raw sample ranges:")
# for begin, end in ranges:
#    print(f"  [{begin}, {end})")

unique_ranges = sorted(set(ranges), key=lambda x: (x[0], x[1]))

print("Sample ranges (sorted):")
print(" idx | begin    | end      | span")
print("-----+----------+----------+---------")
for i, (begin, end) in enumerate(unique_ranges):
   span = end - begin
   print(f" {i:>3} | {begin:>8} | {end:>8} | {span:>7}")

if unique_ranges:
   min_begin = unique_ranges[0][0]
   max_end = max(r[1] for r in unique_ranges)
   overlaps_or_backwards = 0
   gaps = 0
   for (b0, e0), (b1, e1) in zip(unique_ranges, unique_ranges[1:]):
      if b1 < e0:
         overlaps_or_backwards += 1
      elif b1 > e0:
         gaps += 1

   print()
   print(f"Total unique ranges: {len(unique_ranges)}")
   print(f"Global covered window: [{min_begin}, {max_end})")
   print(f"Continuity checks: overlaps/backwards={overlaps_or_backwards}, gaps={gaps}")