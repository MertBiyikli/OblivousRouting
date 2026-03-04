--- a/run_experiments.sh
+++ b/run_experiments.sh
@@ -1,11 +1,15 @@
 #!/usr/bin/env bash
 set -euo pipefail
 
 OUT_CSV="results/results.csv"
 OUT_DIR="${OUT_DIR:-results/logs}"
 TIMEOUT="${TIMEOUT:-600m}"
 
 BIN="${BIN:-./build/oblivious_routing}"
 
+# Optional: propagate your local OR-Tools / Boost / Protobuf runtime libs
+# (harmless if not needed; helps if ldd shows "not found")
+export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
+
 SOLVERS=""
 DATASET=""
 
@@ -96,6 +100,13 @@
 mkdir -p "$(dirname "$OUT_CSV")"
 mkdir -p "$OUT_DIR"
 
+# Sanity: binary exists + executable
+if [[ ! -x "$BIN" ]]; then
+  echo "ERROR: binary not found or not executable: $BIN"
+  echo "Hint: build first: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
+  exit 1
+fi
+
 CSV="$OUT_CSV"
 
 if [[ ! -f "$CSV" ]]; then
@@ -169,26 +180,25 @@
       echo "[RUN] $rel_path | solver=$solver | demand=$demand | timeout=$TIMEOUT"
 
-      # Build argv for the container entrypoint:
-      # - If demand was NOT provided, pass only: <solver> <graph>
-      # - If demand WAS provided, pass: <solver> <graph> <demand>
-      cmd=( "$solver" "graphs/$rel_path" )
+      # Build argv for the native binary:
+      # - If demand was NOT provided, pass only: <solver> <graph>
+      # - If demand WAS provided, pass: <solver> <graph> <demand>
+      cmd=( "$BIN" "$solver" "$g_abs" )
       if [[ "$DEMAND_PROVIDED" -eq 1 ]]; then
         cmd+=( "$demand" )
       fi
 
       status="OK"
       if "$TIMEOUT_BIN" --signal=SIGTERM --kill-after=30s "$TIMEOUT" \
-        docker run --rm --runtime=runc \
-          -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
-          -v "$DATASET_ROOT":/app/graphs \
-          "$IMAGE" \
-          "${cmd[@]}" \
-        > "$log" 2>&1
+        env OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
+          "${cmd[@]}" \
+          > "$log" 2>&1
       then
         status="OK"
       else
         ec=$?
         if [[ "$ec" -eq 124 || "$ec" -eq 137 || "$ec" -eq 143 ]]; then
           status="TIMEOUT"
         elif [[ "$ec" -eq 139 ]]; then
           status="SEGFAULT"
         else
           status="ERROR_$ec"
         fi
       fi