
import pytest

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