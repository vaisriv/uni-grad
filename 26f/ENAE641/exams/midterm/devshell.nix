{ pkgs, perSystem, ... }:
perSystem.devshell.mkShell {
    name = "enae641-midterm formula-sheet devshell";

    packages = with pkgs; [
        # WARN: python3 is overriden to unstable
        (perSystem.nixpkgs-unstable.python3.withPackages (
            ps: with ps; [
                uv
            ]
        ))

        tinymist
    ];

    devshell.startup = {
        activate-venv.text = "source .venv/bin/activate";
    };
}
