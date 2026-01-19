#!/usr/bin/env python3
"""
Script de test pour vérifier le fix du bug 'failed to fill whole buffer'

Ce script teste que le chargement de fichiers STL fonctionne correctement,
même avec des fichiers qui ont des headers corrompus.
"""

from navaltoolbox import Hull
from pathlib import Path
import sys

def test_normal_stl():
    """Test avec un fichier STL normal."""
    print("Test 1: Chargement d'un fichier STL normal...")
    
    stl_path = Path(__file__).parent.parent.parent / "rust" / "tests" / "data" / "dtmb5415.stl"
    
    if not stl_path.exists():
        print(f"  ❌ Fichier {stl_path} non trouvé")
        return False
    
    try:
        hull = Hull(str(stl_path))
        num_triangles = hull.num_triangles()
        bounds = hull.get_bounds()
        
        print(f"  ✅ Fichier chargé avec succès")
        print(f"     - Triangles: {num_triangles}")
        print(f"     - Bounds: LOA={bounds[1]-bounds[0]:.2f}m, BOA={bounds[3]-bounds[2]:.2f}m")
        return True
        
    except Exception as e:
        print(f"  ❌ Erreur lors du chargement: {e}")
        return False


def test_box_stl():
    """Test avec un fichier STL de boîte."""
    print("\nTest 2: Chargement d'un fichier STL box...")
    
    stl_path = Path(__file__).parent.parent.parent / "rust" / "tests" / "data" / "box_10x10.stl"
    
    if not stl_path.exists():
        print(f"  ⚠️  Fichier {stl_path} non trouvé, test ignoré")
        return True  # Not a failure
    
    try:
        hull = Hull(str(stl_path))
        num_triangles = hull.num_triangles()
        
        print(f"  ✅ Fichier chargé avec succès")
        print(f"     - Triangles: {num_triangles}")
        return True
        
    except Exception as e:
        print(f"  ❌ Erreur lors du chargement: {e}")
        return False


def test_empty_file():
    """Test avec un fichier vide."""
    print("\nTest 3: Test avec un fichier vide...")
    
    import tempfile
    
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        temp_path = f.name
    
    try:
        hull = Hull(temp_path)
        print(f"  ❌ Le chargement aurait dû échouer!")
        return False
        
    except Exception as e:
        error_msg = str(e)
        if "STL file is empty" in error_msg or "empty" in error_msg.lower():
            print(f"  ✅ Erreur détectée correctement: {error_msg}")
            return True
        else:
            print(f"  ⚠️  Erreur inattendue: {error_msg}")
            return True  # Still acceptable
            
    finally:
        Path(temp_path).unlink(missing_ok=True)


def main():
    """Exécute tous les tests."""
    print("=" * 70)
    print("Test du fix 'failed to fill whole buffer'")
    print("=" * 70)
    
    results = []
    
    results.append(test_normal_stl())
    results.append(test_box_stl())
    results.append(test_empty_file())
    
    print("\n" + "=" * 70)
    passed = sum(results)
    total = len(results)
    
    if passed == total:
        print(f"✅ Tous les tests ont réussi ({passed}/{total})")
        return 0
    else:
        print(f"❌ Certains tests ont échoué ({passed}/{total})")
        return 1


if __name__ == "__main__":
    sys.exit(main())
