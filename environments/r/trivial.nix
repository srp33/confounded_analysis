with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "trivial-test";
  builder = "/bin/sh";
  args = [ "-c" "echo Hello, Nix! > $out" ];
}
