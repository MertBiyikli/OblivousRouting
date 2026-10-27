
import pytest
from concurrent.futures import ThreadPoolExecutor
from e_routing.scenarios.repository import ScenarioRepository


def test_save_and_retrieve(tmp_path):
    repository = ScenarioRepository(tmp_path)
    result = {"scenario_id": "test-123", "status": "completed"}

    repository.save(result)

    assert repository.get("test-123") == result


def test_missing_scenario(tmp_path):
    repository = ScenarioRepository(tmp_path)

    assert repository.get("nonexistent") is None


def test_reject_duplicate_id(tmp_path):
    repository = ScenarioRepository(tmp_path)
    result = {"scenario_id": "test-123"}

    repository.save(result)

    with pytest.raises(FileExistsError):
        repository.save(result)


def test_reject_invalid_id(tmp_path):
    repository = ScenarioRepository(tmp_path)

    with pytest.raises(ValueError):
        repository.get("../secret")




def test_reserve_and_complete(tmp_path):
    repository = ScenarioRepository(tmp_path)
    repository.reserve("scenario-1")

    assert repository.get("scenario-1")["status"] == "running"

    result = {"scenario_id": "scenario-1", "status": "completed", "metrics": {"congestion": 2.5}}
    repository.complete(result)

    assert repository.get("scenario-1") == result


def test_duplicate_reservation(tmp_path):
    repository = ScenarioRepository(tmp_path)
    repository.reserve("scenario-1")

    with pytest.raises(FileExistsError):
        repository.reserve("scenario-1")




def test_concurrent_reservation(tmp_path):
    repository = ScenarioRepository(tmp_path)

    def reserve():
        try:
            repository.reserve("shared-id")
            return True
        except FileExistsError:
            return False

    with ThreadPoolExecutor(max_workers=10) as pool:
        outcomes = list(pool.map(lambda _: reserve(), range(10)))

    assert outcomes.count(True) == 1
    assert outcomes.count(False) == 9
    assert repository.get("shared-id")["status"] == "running"


def test_failed_serialization_preserves_reservation(tmp_path):
    repository = ScenarioRepository(tmp_path)
    repository.reserve("scenario-1")

    invalid_result = {
        "scenario_id": "scenario-1",
        "status": "completed",
        "metrics": {"congestion": float("nan")},
    }

    with pytest.raises(ValueError):
        repository.complete(invalid_result)

    assert repository.get("scenario-1") == {
        "scenario_id": "scenario-1",
        "status": "running",
    }

    assert list(tmp_path.glob("*.tmp")) == []